#!
rm -rf builddir
meson setup builddir 
meson compile -C builddir/
meson install  -C builddir/