# ChemBlender-Plugin
ChemBlender is a Blender-based plugin for molecular structure visualizaiton that supports 3D modeling of organic molecules, crystal structures, and biomolecules.
To use this plugin, you need to install additional Python packages into the folder X:/XXX/Blender Foundation/Blender 4.X/4.X/python/lib/site-packages, including rdkit, gemmi and biotite.
You can download these ex-packages online or use Blender's Console directly. Here is the method to install Python packages in Blender Console, using rdkit as an example.
>>> import sys
>>> import subprocess
>>> python_exe = sys.executable
>>> subprocess.call([python_exe, '-m', 'ensurepip'])
>>> subprocess.call([python_exe, '-m', 'pip', 'install', '--upgrade', 'pip'])
>>> subprocess.call([python_exe, '-m', 'pip', 'install', 'rdkit'])
This may take a few minutes.
