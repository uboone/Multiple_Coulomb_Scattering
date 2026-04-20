mcs_example: mcs_example.cxx
	`root-config --cxx --cflags` -Isrc -o mcs_example mcs_example.cxx `root-config --libs` -Lsrc
