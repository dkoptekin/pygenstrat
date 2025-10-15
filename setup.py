from setuptools import setup

setup(
    name="pygenstrat",
    version="1.0",
    packages=["pygenstrat"],
    package_dir={"pygenstrat": "src"},
    install_requires=[
        "numpy>=2.0.2",
        "pandas>=2.2.3",
        "psutil>=7.0.0",
        "python-dateutil>=2.9.0",
        "tqdm>=4.67.1",
        "typing_extensions>=4.12.2",
        "scipy>=1.13.1",
        "intervaltree>=3.1.0"
    ],
    setup_requires=[
        "setuptools>=76.0.0",
    ],
    entry_points={
        "console_scripts": [
            "pygenstrat=pygenstrat.main:main",
        ],
    },
    author="Dilek Koptekin",
    author_email="dilek.koptekin@unil.ch",
    description="EIGENSTRAT tools optimized for large datasets",
    url="https://github.com/dkoptekin/pygenstrat",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
