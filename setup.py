from setuptools import setup, find_packages

setup(
      name='MetaSV',
      version='0.1-alpha',
      description='MetaSV: An accurate and integrative structural-variant caller for next generation sequencing',
      author='Bina Technologies',
      author_email='rd@bina.com',
      url='https://github.com/bioinform/metasv',
      packages = find_packages(),
      install_requires = ["pysam", "pybedtools", "pyvcf"],
      package_data = {"metasv": ["resources/*"]},
      scripts=['metasv.py']
     )
