#ifndef FITS_READER_H
#define FITS_READER_H

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include <godot_cpp/variant/string.hpp>
#include <fitsio.h>
extern "C" {
#include <wcs.h>
}

namespace godot {

class FITSReader : public RefCounted {
    GDCLASS(FITSReader, RefCounted);

private:
    fitsfile* fptr;      // FITS file pointer
    int status;          // Status for CFITSIO operations
    int bitpix;          // Bits per pixel
    int naxis;           // Number of axes
    long naxes[2];       // Size of each axis
    bool loaded;         // Flag indicating if file is loaded
    wcsprm* wcs_ptr;     // WCS (World Coordinate System) pointer    
    int nreject;
    // WCS (World Coordinate System) parameters
    double crpix1;       // Reference pixel X
    double crpix2;       // Reference pixel Y
    double crval1;       // Reference value X
    double crval2;       // Reference value Y
    double cd1_1;        // Transformation matrix elements
    double cd1_2;
    double cd2_1;
    double cd2_2;

protected:
    static void _bind_methods();

public:
    FITSReader();
    ~FITSReader();

    bool load_fits(const String& path);
    PackedFloat32Array get_image_data();
    Dictionary get_header_info();
    Dictionary get_info();
    Dictionary get_table_data(int hdu_index);
    Vector2 world_to_pixel(double ra, double dec);
    Vector2 pixel_to_world(double x, double y);
    Dictionary read_spectrum(const String &filename, const String &column_name);
};

}  // namespace godot

#endif // FITS_READER_H

