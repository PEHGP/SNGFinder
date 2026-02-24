from setuptools import setup

setup(
    name="SNGFinder",
    version="0.2.0",
    description="Identifying new genes based on the syntenic method.",
    author="likuan",
    author_email="396777306@qq.com",
    py_modules=["SNGFinder", "SupplementSrc"],
    entry_points={
        "console_scripts": [
            "SNGFinder=SNGFinder:main",
        ],
    },
    python_requires=">=3.12",
)
