from setuptools import setup, find_packages

setup(
    name='dusti',
    version='0.1',
    description='Duplication aware Species Tree Inference',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Megan L. Smith',
    author_email='ms4438@msstate.edu',
    url='https://github.com/SmithLabBio/dusti',
    packages=find_packages(),
    install_requires=[
        'numpy', 'argparse'
    ],
    entry_points={
        'console_scripts': [
            'dusti = dusti.inference:main',  # Command to run from CLI
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)
