from setuptools import setup, find_packages

setup(
    name="cellmerger",                     
    version="0.1.0",                       
    author="Jag Balan",                    
    author_email="balan.jagadheshwar@mayo.edu", 
    description="A short description of your package",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/my_package",  # Replace with your GitHub repo
    packages=find_packages(),              # Automatically find package directories
    install_requires=[
        # List your dependencies here, e.g., 'numpy>=1.18.0'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Replace with your license
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',               # Specify Python version
    entry_points={
        'console_scripts': [
            'my_script=my_package.script:main',  # Command-line entry point
        ],
    },
)

