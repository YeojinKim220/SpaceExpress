# from distutils.core import setup
# from setuptools import find_packages

# # hellow

# with open('requirements.txt') as f:
#     required = f.read().splitlines()

# setup(
#   name = 'spaceepxress',
#   packages=find_packages(include=('SpaceExpress/*',)),
#     package_dir={
#         'SpaceExpress': 'SpaceExpress',
#     },
#   version = '1.0.0',      
#   license='MIT',        
#   description = '',   
#   author = 'Yeojin Kim',                   
#   author_email = 'ykim3030@gatech.edu',      
#   url = 'https://github.com/yeojinkim220',   
#   keywords = ['SpaceExpress', 'Python'],   
#   install_requires=required,
#   classifiers=[
#     # 'Development Status :: 3 - Alpha',      
#     # 'Intended Audience :: Developers',      
#     # 'Topic :: Software Development :: Build Tools',
#     # 'License :: OSI Approved :: MIT License',   
#     # 'Programming Language :: Python :: 3',      
#     # 'Programming Language :: Python :: 3.4',
#     # 'Programming Language :: Python :: 3.5',
#     # 'Programming Language :: Python :: 3.6',
#   ],
#   long_description=open('README.md').read(),
#   long_description_content_type='text/markdown',
# )
from setuptools import setup, find_packages

setup(
    name='SpaceExpress',
    version='0.1.0',
    packages=find_packages(where='src'),  # Automatically find packages in 'src'
    package_dir={'': 'src'},  # Specify the root package directory
    author='Yeojin Kim',
    author_email='ykim3030@gatech.edu',
    url='https://github.com/yeojinkim220',
    keywords=['SpaceExpress', 'Python'],
)