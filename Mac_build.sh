

install_name_tool -change  /opt/X11/lib/libpng16.16.dylib @executable_path/libpng16.16.dylib SEGMENT; 
install_name_tool -change  /Users/tuannt8/Downloads/armadillo-5.000.1/libarmadillo.5.dylib @executable_path/libarmadillo.5.dylib SEGMENT; 
install_name_tool -change  /opt/X11/lib/libX11.6.dylib @executable_path/libX11.6.dylib SEGMENT; 
install_name_tool -change  /Library/Frameworks/GEL.framework/Versions/A/GEL @executable_path/GEL SEGMENT; 