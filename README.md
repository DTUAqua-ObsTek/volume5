A Monte Carlo method can be used to determine the volume of complex or irregularly shaped objects or regions in situations where analytical solutions are not readily available or feasible. This method is widely used in various scientific fields, such as physics, engineering, and computer graphics, where accurate volume determination is required for simulations, modeling, or design purposes. In summary, the object of interest is enclosed within a bounding shape, i.e., in our case a cube, then random points are uniformly distributed within the bounding shape. The number of points falling inside the object is compared to the total number of points generated for the cube. The ratio of the volume of the object to the volume of the cube is proportional to the ratio of points falling inside the object to the total number of points. Thus, by multiplying this ratio by the volume of the bounding shape, an estimate of the object's volume can be obtained. The more points generated, the more accurate the volume estimation becomes. However, even with a relatively small number of points, the Monte Carlo method can provide a reasonable approximation, especially for complex or irregular geometries that are challenging to analyze using traditional mathematical approaches.

The code implements such method for estimation of volume given observed tracks of marine organisms (hence 3D).

The original version of the code is described in:

```
@article{bianco2014analysis,
  title={Analysis of self-overlap reveals trade-offs in plankton swimming trajectories},
  author={Bianco, Giuseppe and Mariani, Patrizio and Visser, Andre W and Mazzocchi, Maria Grazia and Pigolotti, Simone},
  journal={Journal of The Royal Society Interface},
  volume={11},
  number={96},
  pages={20140164},
  year={2014},
  publisher={The Royal Society}
}
```

# Setup

Clone the code repository.

`git clone https://github.com/DTUAqua-ObsTek/volume5.git && cd volume5/`

Compile with the math standard libray:

`gcc volume5.c -o vol5 -lm`

# Usage

`./vol5 2000000000 track.txt > out.out`

Input: Track.txt contains x, y, z values of the track.

Output: out.out contains the radius considered, the volume calculated, the ratio between volume and total volume of the box.

