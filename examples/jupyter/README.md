# Interactive GLVis Jupyter Widget

<!-- Badges generated at https://mybinder.readthedocs.io/en/latest/howto/badges.html -->
[![badge](examples/basic.svg "Basic GLVis + Jupyter Example")](https://mybinder.org/v2/gh/GLVis/pyglvis/HEAD?filepath=examples%2Fbasic.ipynb)
[![badge](examples/plot.svg "Plot grid functions")](https://mybinder.org/v2/gh/GLVis/pyglvis/HEAD?filepath=examples%2Fplot.ipynb)
[![badge](examples/ex1.svg "MFEM's Example 1")](https://mybinder.org/v2/gh/GLVis/pyglvis/HEAD?filepath=examples%2Fex1.ipynb)
[![badge](examples/ex9.svg "MFEM's Example 9")](https://mybinder.org/v2/gh/GLVis/pyglvis/HEAD?filepath=examples%2Fex9.ipynb)

This repository contains a [Jupyter](https://jupyter.org/) widget for the [GLVis](https://glvis.org/) finite element
visualization tool based on the [glvis-js](https://github.com/GLVis/glvis-js) JavaScript/WebAssembly library.

## Usage

```python
from glvis import glvis

glvis(data[, width=640, height=480])

# or assign if you want to update later
g = glvis(data)
# run a cell with `g` to show it
g
```

The `data` object and be one of:

- a `str`, in the format of `*.saved` files
- a `Mesh`, defined in [PyMFEM](https://github.com/mfem/pymfem)
- a `(Mesh, GridFunction)` tuple, defined in [PyMFEM](https://github.com/mfem/pymfem)

[PyMFEM](https://github.com/mfem/pymfem) can be installed with `pip install mfem --no-binary mfem`.


Once you have a `glvis` object there are a few methods that can used to update the
visualization:
```python
# show a new Mesh/GridFunction, resets keys
g.display(data)
# show an updated visualization with the same `Mesh` and `GridFunction`
# dimensions, preserves keys
g.update(data)
# change the image size
g.set_size(width, height)
# force the widget to render. if the widget isn't the last statement in a cell it
# will not be shown without this. see ex9.ipynb
g.show()
```

See the [examples](examples/) directory for additional examples. To test those locally, start a Jupyter notebook server with

```
jupyter notebook
```

## Installation

The GLVis Jupyter widget can be simply installed with `pip`:

```
pip install glvis
```

It order for the installation to be useful you must enable the extension for one or both
of the [Classic Notebook](https://jupyter-notebook.readthedocs.io/en/stable/) and
[Jupyter Lab](https://jupyterlab.readthedocs.io/en/stable/), see the next two sections:

### Jupyter Notebook

To use the widget with the basic Notebook enable it with `jupyter nbextension enable`:

```
jupyter nbextension enable --py glvis
```

After enabling the extension you can verify you're good to go with:

```
jupyter nbextension list
```

The output should be something like:

```
Known nbextensions:
  config dir: path/to/nbconfig
    notebook section
      glvis-jupyter/extension  enabled
      - Validating: OK
  <possibly a different config dir>
      jupyter-js-widgets/extension  enabled
      - Validating: OK
```

If `glvis-jupyter` and `jupyter-js-widgets` are not both listed, try the following:

```
jupyter nbextension install --user --py glvis
jupyter nbextension enable --user --py glvis
jupyter nbextension install --user --py widgetsnbextension
jupyter nbextension enable --user --py widgetsnbextension
```
You may also need to run these if you have upgraded to a newer version of the GLVis Jupyter widget.

### Jupyter Lab

[JupyterLab](https://jupyterlab.readthedocs.io) requires another set of install commands:

```
jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
jupyter labextension install glvis-jupyter
```

## Development

Development installation requires `npm`, after installing:

```
git clone https://github.com/glvis/pyglvis.git
cd pyglvis
pip install -e .
```


### Developing in Jupyter Notebook

```
jupyter nbextension install --py --symlink --sys-prefix glvis
jupyter nbextension enable --py --sys-prefix glvis
```

### Developing in Jupyter Lab

```bash
jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
# I believe you need node in the path Lab uses for this to work, I see an extension load error
# in a context where I don't have it:
# Failed to load resource: the server responded with a status of 500 (Internal Server Error)
#   lab/api/extensions?1610138347763
# Which is just a python stacktrace, ending with:
#   raise ValueError(msg)
#   ValueError: Please install Node.js and npm before continuing installation.
jupyter labextension link ./js
```


### Troubleshooting

If you run into errors related to node/npm that aren't helpful try:

```bash
cd pyglvis
make clean
cd js
# fix errors in these steps, run `make -C .. clean` each time
npm install
npx webpack
```

## Releasing

### Releasing a new version of glvis-jupyter on NPM:

- Update the required version of `glvis` in `js/package.json`

- Update the version in `js/package.json`

```bash
# clean out the `dist` and `node_modules` directories
git clean -fdx
npm install
npm publish
```

### Releasing a new version of glvis on PyPI:

- Update `glvis/_version.py`
   - Set release version
   - Update `extension_version` to match `js/package.json`

- `git add` and `git commit` changes
  - `glvis/_version.py`, `js/package.json`, and `js/package-lock.js`


You will need [twine](https://pypi.org/project/twine/) to publish to PyPI, install with `pip`.

```
python setup.py sdist bdist_wheel
twine upload dist/*
git tag -a X.X.X -m 'comment'
git push --tags
```
