'ExtracSolenB_BminBecr0.xx.csv' contain the axial magnetic field from the
etraction zperture (z=0mm) to z=400mm. The axial field extent may need to
be extended further to fully capture the canonical momentum contribution
to beam emittance due to the source extraction solenoidal fields.


'./Electrode/SuSI_ExtracAperxxmm' contain the electrode geometries for 
extraction aperture of diameter of d=xx mm. The file'Extraction.dat' 
contains:

>>> electrode_data
{'cond_ids': [1, 2, 3], 'voltages': [24000.0, -2000.0, 0.0], 'cond_names': ['Source', 'Puller', 'Ground']}

Whereas the actual geometry of the extraction electrodes are contained in the .wob file.


Work-in-progress
================
- push the files and instructions necessary to generate axial field
profile and electrode geometry files

