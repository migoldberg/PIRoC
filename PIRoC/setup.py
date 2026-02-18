from __future__ import absolute_import
from pathlib import Path
from setuptools import find_packages, setup

ROOT = Path(__file__).resolve().parent

readme_path = ROOT / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="PIRoC",
    version="0.1.0",
    author="Mark Goldberg",
    license="MIT",
    description=(
        "A tool for decontamination of inverterbrate transcriptomes using phylogenetic trees"
    ),
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={"console_scripts": ["PIRoC = PIRoC.__main__:main"]},
    install_requires=["ete3>=3.1.0"],
    python_requires=">=3.8",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/markgoldberg11/dev-PIRoC",
    keywords=[
        "bioinformatics",
        "phylogenetics",
        "contamination",
        "orthogroups",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
