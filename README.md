# Microfacet Materials Renderer

This project is based on a class project of CS 190I Introduction to Offline Rendering at UC Santa Barbara. In this project, I
implemented the BRDF evaluation function with the Fresnel term for reflection on metallic objects. I also computed the 
Beckmann normal distribution function, which is similar to a Gaussian distribution, that defines how the 
microfacets' normals are distributed. I also adopted the importance sampling microfacet BRDF that works better
with Beckmann NDF.

## Usage
After you have cloned or downloaded the code, enter the folder/repo and you can use the following to compile and
run the code:
```
$ mkdir build
$ cd ./build
$ cmake ..
$ make
```
After you have successfully compiled the code, you can use `./pathtracer` with appropriate arguments to execute
the code. For example, the command below will produce an image **out.png** of size 500 * 380 with 64 samples per
pixel using the **CBdragon_au.dae** scene.

`./pathtracer -r 500 380 -f out.png -s 64 ../scenes/CBdragon_au.dae`

This repo comes with two scenes: **CBbunny_cu.dae** and **CBdragon_au.dae**.

## Sample Output
Here're the generated results with 128 ssp for both the scenes:

![alt text](https://github.com/gitzhujiarui/microfacet-materials-renderer/blob/main/output/bunny.png?raw=true)

![alt text](https://github.com/gitzhujiarui/microfacet-materials-renderer/blob/main/output/dragon.png?raw=true)
