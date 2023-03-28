from setuptools import setup
# NOTE: necessary to import __version__

from os import path
with open(path.join(path.abspath(path.dirname(__file__)), "README.md")) as f:
    readme = f.read()
setup(
    name="dcca",
    version="0.2.0",
    description="Python implementation of the Time-Lagged Detrended Cross-Correlation Coefficient Analysis (DCCA)",
    url="https://github.com/LeonardoAlchieri/dcca",
    author="Leonardo Alchieri",
    author_email="leonardo@alchieri.eu",
    license="GPLv3",
    packages=["dcca"],
    install_requires=["numpy"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        # 'Programming Language :: Python :: 3.4',
        # 'Programming Language :: Python :: 3.5',
        # 'Programming Language :: Python :: 3.6',
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
