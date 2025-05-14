from setuptools import setup, find_packages

setup(
    name="deeplotyper",
    version="0.1.0",
    description="Genomic ↔ transcript coordinate mapping and haplotype remapping tools",
    author="Eli Niktab",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.78"
    ],
    python_requires=">=3.7",
)