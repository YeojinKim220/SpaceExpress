from setuptools import setup, find_packages

setup(
    name='spaceexpress',
    packages=find_packages(where='src'),  # Automatically find packages in 'src'
    package_dir={'': 'src'},  # Specify the root package directory
    author='Yeojin Kim',
    author_email='ykim3030@gatech.edu',
)
