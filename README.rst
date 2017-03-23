
Neuropype_graph: Neuropype project for graph analysis, can be used from neuropype_ephy and nipype

Warning, should add all the directories from radatools to the PATH env variable

1) Download radatools sotware:

http://deim.urv.cat/~sergio.gomez/radatools.php#download

2) Download and extract the zip file

3) Add following lines in your .bashrc:

RADA_PATH=/home/david/Tools/Software/radatools-4.0-linux64

(replace /home/david/Tools/Software by yout path to radatools)

export PATH=$PATH:$RADA_PATH/Network_Tools
export PATH=$PATH:$RADA_PATH/Network_Properties
export PATH=$PATH:$RADA_PATH/Communities_Detection 
export PATH=$PATH:$RADA_PATH/Communities_Tools   
