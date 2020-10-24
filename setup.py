import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = {}
with open("iqkm/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    name="iqkm",
    version=version["__version__"],
    author="Jingdi Li",
    author_email="lijingdioo@outlook.com",
    description="Identification and quantification of KEGG Modules in metagenomes/genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Finn-Lab/iqKM/",
    entry_points={"console_scripts": ["iqkm = iqkm.__main__:main"]},
    install_requires=["networkx", "pysam", "markdown", "pytest", "pylint"],
#    dependency_links=['https://github.com/snayfach/MicrobeCensus'],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3.8.3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
    ],
)

