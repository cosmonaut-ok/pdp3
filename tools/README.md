# Tools

### Quick description

- **image_2_et_component.py** - Build image (png) of two 2D plots of electrostatical potential, related to time (E_r/E_z, t coordinates)
- **image_2e_zt_component.py** - Build image (png) of two 2D color maps of electrostatical potential, dependency from time and z-coordinate (z, t coordinates, E_r/E_z - color maps)
- **images_3_component.py** - Build image (png) of three color maps of electrostatical potential E_r/E_z and particle beam density \rho_{beam}, dependent on r and z coordinates in selected time moment
- **movie_3_component.py** - Build movie of three color maps of electrostatical potential E_r/E_z and particle beam density \rho_{beam}, dependent on r and z. Shows, how E_r, E_z and \rho_{beam} changed in time.
- **quick_parameters_calculator.py** - quickly calculate some important parameters of model (accepts parameters.xml as an argument)
- **view_2_et_component.py** - same as image_2_et_component.py, but outputs plot on screen interactively, instead of writing to file(s)
- **view_2e_zt_component.py** - same as image_2e_zt_component.py, but outputs plot on screen interactively, instead of writing to file(s)
- **view_3_component.py** - same as images_3_component.py and movie_3_component.py, but outputs plot on screen interactively, instead of writing to file(s)

### Usage

All tools uses pdp3 configuration file `parameters.xml` as the source of truth. So, all of the tools accepts required argument `properties_path`, which is path to `parameters.xml` file. Other options (depends on tool):

- **-h/--help** - __view help for any tool__

- **--timestamp** - set moment of time which plot should be built for
- **--time-range** - set time range which images with plot, movie, or interactive animation should be generated for
- **--data-set-range** - same as time range, but first and last data sets gives directly
- **--cmap** - color map, which should be used for color map plots ( [cmap list](https://matplotlib.org/examples/color/colormaps_reference.html) )
- **--clim-e-r/--clim-e-z (etc)** - set color limits for E_r, E_z etc *(calculated automatically by default)*
- **--beam-scale-factor** - color limits for particle beam calculated automatically. User can set scale factor to highlight some details.
- **radius** - set radius of point for feld parameters measurement *(used in image_2_et_component.py etc.)*
- **longitude** - set longitude of point for feld parameters measurement *(used in image_2_et_component.py etc.)*
- **--view** - view plots interactively, while processing movie or images set *(used in image_.. and movie_.. tools)*
- **--dry-run** - do not save any files. Just simulate *(used for debug in image_.. and movie_.. tools)*
- **--video-file** - set path to video file *(used in movie_.. tools)*
- **--images-path** - set path to images set *(used in images_3.. tools)*
- **--image-file** - set path to image file *(used in image_2.. tools)*

### Examples

``` shell
$ ./tools/images_3_component.py /home/user/pdp3_model/parameters.xml --cmap=terrain --clim-e-r=0:5e5 --clim-e-r=0:5e5
```

``` shell
./tools/image_2_et_component /home/user/pdp3_model/parameters.xml --radius=0.01 --longitude=0.1 --time-range=1e-9:2.7e-9
```
