# warp_ecr
Warp slice simulations of ECR ion sources and extraction beamlines at MSU using highly reduced models.

Code scripts assembled by Alfonse Pham based reduced ECR ion source models by 
Damon Todd and Daniela Leitner (working at LBNL) which were extended by Daniel Winklehner (while a grad student at MSU). 
[Alfonse ... please extend description here]

Contact Info for Present Collaborators:

Alfonse Pham 
Michgan State University 
National Superconducting Cyclotron Laboratory
pham@nscl.msu.edu

Daniel Winklehner 
[Alfonse put info]

Derek Neben, Chun Yan (Jonathan) Wong, Daniela Leitner 
[Alfonse put info ... here is an edit] 

Steven M. Lund
Physics and Astronomy Department
Facility for Rare Isotope Beams
Michigan State University
lund@frib.msu.edu
517-908-7291

To initialize the repository, 

   % git clone https://github.com/smlund/warp_ecr

This will create a directory, ./warp_ecr where command was 
run with the archive files.   

To get the latest version descend into directory warp_ecr and run:

  % git pull 

When modifying the repository (for those with edit privilege, please contact 
me if you want to contribute and I will add you) 

  ... edit files etc then checkin (use README.md as example here) using: 
  % git add README.md 
  % git commit -m "SML: updated README.md file" 
  % git push

You may be prompted for user name and password info in this step.  

To remove a file "file" from git control (to not include in future pulls), use 

  % git rm file
  % git commit -m "SML: removing file from repo"
  % git push 

The local copy of the file removed can be retained on disk by subsituting

  % git rm --cached file 

in the above. Another way to remove a file from the master node repo is to
go to the github web side, click on a file, click the delete button and then
confirm at the bottom of the page.   

[Alfonse ... please update below as we check in needed field description files]
Field files needed for the lattices described by the simulation are stored in the 
Dropbox file sharing system under an account of Steven Lund. Code users do NOT need 
dropbox access to retrieve the needed field description files for lattice elements.
A unix shell script "ecr-lat-fields-dropbox-fetch" is provided that
employs wget with dropbox links to download the needed field element
description files in various lat_element_name subdirectories.  Examples:
     lat_susi        Venus source 
	   lat_s4          s4 solenoid
	   ...
Run this executable script using:
 % ./ecr-lat-fields-dropbox-fetch

This script and the linked Dropbox account must be maintained consistently
as lattice element descriptions change. Code users may want access to the
Dropbox account which contains much info/input on the generation of lattice
element field files, code input (Laplace, Poisson, ...), element plots,
etc. Contact Steve Lund for access if you have a dropbox account and this
account can be "shared" with you. Note that this directories is large due to 
the size of field element arrays and shared account contribute to a dropbox 
users quota. So a larger (paid) dropbox account would likely be needed for 
full sharing. 

