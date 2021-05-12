from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="dbspro",
    author="Tobias Frick, Pontus Höjer",
    url="https://github.com/TobiasFrick/DBS-Pro/",
    description="DBS-Pro pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    package_dir={"": "src"},
    install_requires=[
        "pysam",
        "pandas",
        "dnaio",
        "tqdm",
        "snakemake",
        "importlib_resources; python_version<'3.7'",
        "umi_tools",
    ],
    extras_require={
        "dev": [
            "flake8",
            "pytest",
            ],
    },
    package_data={"dbspro": ["rules.smk", "report_template.ipynb", "dbspro.yaml", "config.schema.yaml",
                             "ABC-sequences.fasta"]},
    packages=find_namespace_packages("src"),
    entry_points={"console_scripts": ["dbspro = dbspro.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
