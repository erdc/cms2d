from conan import ConanFile
from conan.tools.cmake import CMake

class XMDF_Conan(ConanFile):
    name = "xmdf"
    user = "Aquaveo"
    version = "2.2"
    description = "a C and Fortran language library providing a standard format for the geometry data storage of river cross-sections, 2D/3D structured grids, 2D/3D unstructured meshes, geometric paths through space, and associated time data."
    license = "MIT"
    homepage = "https://www.xmswiki.com/wiki/XMDF"
    settings = "build_type"
    export_sources = "src*"
    requires = "zlib/1.3.1", "szip/2.1.1", "hdf5/1.14.5"
    generators = "CMakeDeps", "CMakeToolchain"

    def build(self):
        # Either using some of the Conan built-in helpers
        cmake = CMake(self)
        cmake.configure()  # equivalent to self.run("cmake . <other args>")
        cmake.build() # equivalent to self.run("cmake --build . <other args>")

