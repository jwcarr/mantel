import setuptools

with open("README.md", encoding="utf-8") as file:
    long_description = file.read()

setuptools.setup(
    name="mantel",
    use_scm_version={
        "write_to": "mantel/_version.py",
        "write_to_template": '__version__ = "{version}"',
        "fallback_version": "???",
    },
    author="Jon Carr",
    author_email="jcarr@sissa.it",
    description="Python implementation of the Mantel test, a significance test of the correlation between two distance matrices",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jwcarr/mantel",
    license="MIT",
    packages=["mantel"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.6",
    install_requires=["numpy>=1.10", "scipy>=1.0"],
    setup_requires=["setuptools_scm"],
)
