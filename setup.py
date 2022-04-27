import os
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='trm.utils',
    version='0.0.1',
    packages = ['trm.utils',],

    # metadata
    author='Tom Marsh',
    author_email='t.r.marsh@warwick.ac.uk',

    description="Python utility classes and functions",
    long_description=long_description,

    url='http://www.astro.warwick.ac.uk/',

    # Choose your license
    license='BSD',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Astronomers',
        #'Topic :: Astronomy :: Photometric reduction',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: BSD License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],

    install_requires=['astropy','matplotlib','numpy'],

    # Makes significant use of f-strings which came in python v3.6
    python_requires='>=3.6',

  )
