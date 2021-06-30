#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}


// Diffuse BSDF //
Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return reflectance * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return reflectance * (1.0 / PI);
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1;
  reflect(wo,wi);
  return  reflectance * (1./wo[2]);
}





// Microfacet BSDF ///////////////////////

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    // Shadowing-masking term
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    /*
	  TODO:
	  Compute Beckmann normal distribution function (NDF) here.
	  You will need the roughness "alpha".
	*/
    Vector3D h_normalized = h / h.norm();
    Vector3D normal = Vector3D(0, 0, 1);
    // double cos_theta_h = dot(h_normalized, normal);
    double cos_theta_h = cos_theta(h_normalized);
    if (cos_theta_h == 0)
      return 0;

    // double tan_theta_h = cross(h_normalized, normal).norm() / dot(h_normalized, normal);
    double tan_theta_h = (sqrt(1 - cos_theta_h * cos_theta_h)) / cos_theta_h;
    double result = exp((-1) * tan_theta_h * tan_theta_h / this->alpha / this->alpha) / 
                    (PI * this->alpha * this->alpha * pow(cos_theta_h, 4.0));
	
    return result;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    /*
      TODO:
      Compute Fresnel term for reflection on air-conductor interface.
      You will need both "eta" and "k", both of which are Spectrum.
    */
    Spectrum R_s = ((eta * eta + k * k) - 2 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi)) /
                  ((eta * eta + k * k) + 2 * eta * cos_theta(wi) + cos_theta(wi) * cos_theta(wi));
    Spectrum R_p = ((eta * eta + k * k) * cos_theta(wi) * cos_theta(wi) - 2 * eta * cos_theta(wi) + Spectrum(1, 1, 1)) /
                  ((eta * eta + k * k) * cos_theta(wi) * cos_theta(wi) + 2 * eta * cos_theta(wi) + Spectrum(1, 1, 1));
    Spectrum result = (R_s + R_p) / 2;
    return result;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    /*
      TODO:
      Implement microfacet model here
      Note that you will return the BRDF only, without the cosine term
    */
    Vector3D n = Vector3D(0, 0, 1);

    if (dot(n, wo) <= 0 || dot(n, wi) <= 0)
      return Spectrum();

    Vector3D h = (wo + wi) / 2;
    h.normalize();
    Spectrum result = F(wi) * G(wo, wi) * D(h) / (4 * dot(n, wo) * dot(n, wi));
    
    return result;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    /*
	  *Importance* sample Beckmann normal distribution function (NDF) here.
	  Note: You should fill in the sampled direction *wi and the corresponding *pdf,
	   	    and return the sampled BRDF value.
	*/

    // - Sample theta_h and phi_h
    Vector2D random_nums = sampler.get_sample();
    double r1 = random_nums.x;
    double r2 = random_nums.y;
    double theta_h = atan(sqrt((-1) * alpha * alpha * log(1 - r1)));
    double phi_h = 2 * PI * r2;

    // - Calculate outgoing direction *wi from sampled h
    Vector3D h, outgoing_dir;
    h.x = sin(theta_h) * cos(phi_h);
    h.y = sin(theta_h) * sin(phi_h);
    h.z = cos(theta_h);
    h.normalize();
    // half vector reflector
    *wi = (2.0 * dot(wo, h)) * h - wo;

    // This is left as a gift for you. Make sure to comment this out 
    // when you are using the cosineHemisphereSampler
    if (dot(wo, h) <= 0.0 || wo.z <= 0.0 || wi->z <= 0.0) {
        *pdf = 10.0;
        return Spectrum();
    }

    // - Calculate *pdf of sampling *wi
    // double sin_theta_h = sin_theta(h);
    // double cos_theta_h = cos_theta(h);
    double sin_theta_h = sin(theta_h);
    double cos_theta_h = cos(theta_h);
    double tan_theta_h = cos_theta_h == 0 ? 0 : tan(theta_h);
    double p_theta = 0;
    if (cos_theta_h != 0) {
      p_theta = 2 * sin_theta_h / (alpha * alpha * pow(cos_theta_h, 3.0)) *
                     exp((-1) * tan_theta_h * tan_theta_h / alpha / alpha);
    }
    double p_phi = 1.0 / 2 / PI;
    double p_w_h = sin_theta_h <= 0 ? 0 : (p_theta * p_phi / sin_theta_h);
    double p_w = dot(*wi, h) <= 0 ? 0 : (p_w_h / (4 * dot(*wi, h)));
    if (p_w <= 0)
      p_w = 10.0;
    *pdf = p_w;

    // Comment this line once you have finished with the importance sampling 
    // *wi = cosineHemisphereSampler.get_sample(pdf);
    
    return MicrofacetBSDF::f(wo, *wi);
}

//////////////////////////////////////////




// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  *pdf = 1.0f;

  float cosTheta = cos_theta(wo);
  bool entering = cosTheta > 0;
  float ei = 1.f, et = ior;
  if (!entering) {
      swap(ei, et);
      cosTheta = -cosTheta;
  }
  float inveta = et / ei;
  float inveta2 = inveta * inveta;

  if (refract(wo, wi, ior))
      return inveta2 / cosTheta * transmittance;
  else
      return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // compute Fresnel coefficient and use it as the probability of reflection
  float R0 = (ior-1.0)*(ior-1.0)/((ior + 1.0)*(ior + 1.0));
  float cosTheta = cos_theta(wo);
  float f = 1 - fabs(cosTheta);
  float g = ((f * f) * (f * f)) * f;
  float fresnel_coe = R0 + (1.0 - R0)*g;

  bool entering = cos_theta(wo) > 0;
  float ei = 1.f, et = ior;
  if (!entering) {
      swap(ei, et);
      cosTheta = -cosTheta;  // be careful here, want cosTheta to be
                             // positive for everything below
  }
  float inveta = et / ei;
  float inveta2 = inveta * inveta;

  if (!refract(wo, wi, ior)) {
    // total internal reflection; always reflect
    *pdf = 1.0;
    reflect(wo, wi);
    return (1 / cosTheta) * reflectance;
  }

  if ((double)(std::rand()) / RAND_MAX < fresnel_coe) {
    *pdf = fresnel_coe;
    reflect(wo, wi);
    return (fresnel_coe / cosTheta) * reflectance;
  } else {
    // refraction ray has already been computed
    float one_minus_fresnel = 1.0f - fresnel_coe;
    *pdf = one_minus_fresnel;
    return (one_minus_fresnel * inveta2 / cosTheta) * transmittance;
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0],-wo[1],wo[2]);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  bool entering = cos_theta(wo) > 0;

  float ei = 1.f, et = ior;
  if (!entering) swap(ei, et);

  float sini2 = sin_theta2(wo);
  float eta = ei / et;
  float sint2 = eta * eta * sini2;
  if (sint2 > 1.f) return false;
  float cost = sqrt(1.0f - sint2);

  if (entering) cost = -cost;
  float sint_over_sini = eta;

  *wi = Vector3D(-sint_over_sini * wo.x, -sint_over_sini * wo.y, cost);

  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
