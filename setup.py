import io
import os
import setuptools


def read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_requirements(path):
    content = read(path)
    return [req for req in content.split("\n") if req != "" and not req.startswith("#")]


def get_packages():
    return [
        package
        for package in setuptools.find_packages()
        if package.startswith("comut_plot")
    ]


setuptools.setup(
    name="comut_plot",
    description="Creation of Comutation Plots.",
    version="0.1",
    long_description=read("README.md"),
    author="Philipp Hähnel",
    author_email="phahnel@broadinstitute.org",
    license="MIT",
    url="",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=get_packages(),
    # package_data={'': ['*.gz', '*.txt']},
    include_package_data=True,
    install_requires=get_requirements("requirements.txt"),
)
