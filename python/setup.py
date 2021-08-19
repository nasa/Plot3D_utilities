from setuptools import setup
from Cython.Build import cythonize


_config = {
    "name": "Plot3D",
    "version":"0.1.1",
    "url": "",
    "author": "Paht Juangphanich",
    "author_email": "paht.juangphanich@nasa.gov",    
    "install_requires":['bz2file','tqdm','scipy','pandas','numpy','matplotlib'],
    "packages":["plot3d"],
    'license':"GNU GPLv3",
    'zip_safe':False,
}

def main() -> int:
    """ Execute the setup command.
    """

    def version():
        """ Get the local package version. """
        return _config["version"]

    _config.update({
        "version": version(),
    })

    setup(**_config, ext_modules = cythonize('plot3d/module/plot3d_point_match.pyx',language="3"))
    return 0


# Make the script executable.
if __name__ == "__main__":
    raise SystemExit(main())
