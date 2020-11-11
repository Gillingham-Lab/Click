import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ClickReaction", # Replace with your own username
    version="0.1.1",
    author="Basilius Sauter",
    author_email="basilius.sauter@gmail.com",
    description="A collection of chemical reaction formulations for use with rdkit. Requires rdkit.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Gillingham-Lab/Click",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Utilities",
    ],
    python_requires='>=3.6',
)
