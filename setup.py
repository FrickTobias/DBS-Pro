import sys
from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="DBSpro",
    author="Tobias Frick",
    url="https://github.com/FrickTobias/DBSpro",
    description="DBSpro analysis pipeline",
    long_description_content_type="text/markdown",
    licence="",
    python_requirements=">=3.6",
    package_dir={"":"src"},
    packages=find_packages("src"),
    entry_points={"console_scripts": ["DBSpro = DBSpro.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Licence :: ..."
        "Programming language :: Python :: 3",
        "Programming language :: Python :: 3.6",
        "Programming language :: Python :: 3.7",

    ]
)