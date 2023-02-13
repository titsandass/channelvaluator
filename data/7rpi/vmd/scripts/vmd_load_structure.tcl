#set dir "/home/caverweb/jobs/w5/w5b6cj/caver/input"

mol load pdb ../data/protein.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

