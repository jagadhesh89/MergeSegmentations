from setuptools import setup, find_packages

setup(
    name="cellmerge",                     
    version="0.1.0",                       
    author="Jag Balan",                    
    author_email="balan.jagadheshwar@mayo.edu", 
    description="Merge nuclei segmentations from multiple models",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jagadhesh89/MergeSegmentations/tree/main/cellmerge",  
    packages=find_packages(),              
    install_requires=["tiatoolbox >= 1.4.0","numpy >=1.23.5","annoy >= 1.17.3","matplotlib >= 3.6.2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',               
    entry_points={
        'console_scripts': [
            'cellmerge=cellmerge.cellmerge:main',
        ],
    },
)
