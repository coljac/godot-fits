#include "fits_reader.h"
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/array.hpp>
#include <fitsio.h>
#include <string>
#include <vector>
#include <cstdio>
#include <wcs.h>
#include <wcshdr.h>
// GAH
#include <wcserr.h>

using namespace godot;

void FITSReader::_bind_methods() {
    ClassDB::bind_method(D_METHOD("load_fits", "path"), &FITSReader::load_fits);
    ClassDB::bind_method(D_METHOD("get_info"), &FITSReader::get_info);
    ClassDB::bind_method(D_METHOD("get_image_data"), &FITSReader::get_image_data);
    ClassDB::bind_method(D_METHOD("get_header_info"), &FITSReader::get_header_info);
    ClassDB::bind_method(D_METHOD("get_table_data", "hdu_index"), &FITSReader::get_table_data);
    ClassDB::bind_method(D_METHOD("world_to_pixel", "ra", "dec"), &FITSReader::world_to_pixel);
    ClassDB::bind_method(D_METHOD("pixel_to_world", "x", "y"), &FITSReader::pixel_to_world);
    ClassDB::bind_method(D_METHOD("read_spectrum", "filename", "column_name"), &FITSReader::read_spectrum);
}

FITSReader::FITSReader() : fptr(nullptr), status(0), bitpix(0), naxis(0), loaded(false),
    crpix1(0), crpix2(0), crval1(0), crval2(0),
    cd1_1(1), cd1_2(0), cd2_1(0), cd2_2(1), wcs_ptr(nullptr), nreject(0)
{
    naxes[0] = 0;
    naxes[1] = 0;
}

FITSReader::~FITSReader() {
    if (loaded && fptr) {
        fits_close_file(fptr, &status);
        loaded = false;
    }
    if (wcs_ptr) {
        wcsfree(wcs_ptr);
        delete wcs_ptr;
    }
}

bool initialize_wcsprm(fitsfile *fptr, wcsprm **wcs_ptr) {
    int status = 0, nkeys = 0, nreject = 0, nwcs = 0;
    char *header = nullptr;
    struct wcsprm *wcs_struct = nullptr;

    std::cout << "Reading FITS header..." << std::endl;

    if (fits_hdr2str(fptr, 1, nullptr, 0, &header, &nkeys, &status)) {
        std::cerr << "Error reading FITS header: " << status << std::endl;
        return false;
    }
    std::cout << "FITS Header Dump:" << std::endl;
    std::cout << header << std::endl;

    std::cout << "Header read successfully, number of keys: " << nkeys << std::endl;
    std::cout << "Parsing WCS information using wcspih()..." << std::endl;

    int wcspih_status = wcspih(header, nkeys, 0, 0, &nreject, &nwcs, &wcs_struct);
    free(header);  // Free header immediately after parsing
    if (wcspih_status || nwcs == 0) {
        std::cerr << "Error parsing WCS from FITS header." << std::endl;
        return false;
    }

    std::cout << "WCS parsing completed. Number of WCS structures found: " << nwcs << std::endl;
    std::cout << "Number of rejected keywords: " << nreject << std::endl;
    
    // Assign wcs_struct[0] directly to wcs_ptr, transferring ownership
    *wcs_ptr = wcs_struct;  // Now wcs_ptr points to the allocated structure
// std::cout << "Parsed WCS Parameters Before wcsset():" << std::endl;
// std::cout << "CTYPE: (" << (*wcs_ptr)->ctype[0] << ", " << (*wcs_ptr)->ctype[1] << ")" << std::endl;
// std::cout << "CRPIX: (" << (*wcs_ptr)->crpix[0] << ", " << (*wcs_ptr)->crpix[1] << ")" << std::endl;
// std::cout << "CRVAL: (" << (*wcs_ptr)->crval[0] << ", " << (*wcs_ptr)->crval[1] << ")" << std::endl;
// std::cout << "CD Matrix: [[" << (*wcs_ptr)->cd[0] << ", " << (*wcs_ptr)->cd[1] << "], ["
//                                   << (*wcs_ptr)->cd[2] << ", " << (*wcs_ptr)->cd[3] << "]]" << std::endl;

    // if ((*wcs_ptr)->cd[3] == 0.0) {  // CD2_2 should not be zero
    //     std::cerr << "WARNING: CD2_2 is missing or zero, setting manually!" << std::endl;
    //     (*wcs_ptr)->cd[3] = 0.000277777777777778;  // Set a reasonable default
    //     (*wcs_ptr)->pc[3] = 1.0;  // Set a reasonable default
    // }
// std::cout << "Parsed WCS Parameters Before wcsset():" << std::endl;
// std::cout << "CTYPE: (" << (*wcs_ptr)->ctype[0] << ", " << (*wcs_ptr)->ctype[1] << ")" << std::endl;
// std::cout << "CRPIX: (" << (*wcs_ptr)->crpix[0] << ", " << (*wcs_ptr)->crpix[1] << ")" << std::endl;
// std::cout << "CRVAL: (" << (*wcs_ptr)->crval[0] << ", " << (*wcs_ptr)->crval[1] << ")" << std::endl;
// std::cout << "CD Matrix: [[" << (*wcs_ptr)->cd[0] << ", " << (*wcs_ptr)->cd[1] << "], ["
//                                   << (*wcs_ptr)->cd[2] << ", " << (*wcs_ptr)->cd[3] << "]]" << std::endl;
//     (*wcs_ptr)->flag = -1;
//     // Finalize WCS initialization

std::cout << "Parsed CD Matrix Before wcsset():" << std::endl;
std::cout << "[[" << (*wcs_ptr)->cd[0] << ", " << (*wcs_ptr)->cd[1] << "], ["
                  << (*wcs_ptr)->cd[2] << ", " << (*wcs_ptr)->cd[3] << "]]" << std::endl;


    std::cout << "Parsed WCS Parameters Before wcsset():" << std::endl;
std::cout << "CTYPE: (" << (*wcs_ptr)->ctype[0] << ", " << (*wcs_ptr)->ctype[1] << ")" << std::endl;
std::cout << "CRPIX: (" << (*wcs_ptr)->crpix[0] << ", " << (*wcs_ptr)->crpix[1] << ")" << std::endl;
std::cout << "CRVAL: (" << (*wcs_ptr)->crval[0] << ", " << (*wcs_ptr)->crval[1] << ")" << std::endl;
std::cout << "CD Matrix: [[" << (*wcs_ptr)->cd[0] << ", " << (*wcs_ptr)->cd[1] << "], ["
                                  << (*wcs_ptr)->cd[2] << ", " << (*wcs_ptr)->cd[3] << "]]" << std::endl;
std::cout << "CDELT: (" << (*wcs_ptr)->cdelt[0] << ", " << (*wcs_ptr)->cdelt[1] << ")" << std::endl;
std::cout << "LONPOLE: " << (*wcs_ptr)->lonpole << std::endl;
std::cout << "LATPOLE: " << (*wcs_ptr)->latpole << std::endl;
std::cout << "PC Matrix: [[" << (*wcs_ptr)->pc[0] << ", " << (*wcs_ptr)->pc[1] << "], ["
                                  << (*wcs_ptr)->pc[2] << ", " << (*wcs_ptr)->pc[3] << "]]" << std::endl;
std::cout << "NAXIS: " << (*wcs_ptr)->naxis << std::endl;
// (*wcs_ptr)->cdelt[0] = 0.0;
// (*wcs_ptr)->cdelt[1] = 0.0;
// (*wcs_ptr)->flag = -1;
// if ((*wcs_ptr)->lonpole > 1.0e10 || (*wcs_ptr)->lonpole < -1.0e10) {
    // std::cerr << "WARNING: LONPOLE is an invalid value (" << (*wcs_ptr)->lonpole << "), setting to 180.0" << std::endl;
    // (*wcs_ptr)->lonpole = 180.0;
// }

    wcserr_enable(1);

    int wcsset_status = wcsset(*wcs_ptr);
    std::cout << "wcsset() return status: " << wcsset_status << std::endl;

if (wcsset_status != 0) {
    if ((*wcs_ptr)->err) {
        std::cerr << "WCSLIB ERROR: " << (*wcs_ptr)->err->msg << std::endl;
        std::cerr << "Error occurred in function: " << (*wcs_ptr)->err->function
                  << " at " << (*wcs_ptr)->err->file << ":"
                  << (*wcs_ptr)->err->line_no << std::endl;
    } else {
        std::cerr << "WCSLIB ERROR: wcsset() failed with status " << wcsset_status << " but no error message found." << std::endl;
    }
}




    // Check if WCS was initialized properly
    if ((*wcs_ptr)->flag == 0) {
        std::cerr << "ERROR: WCS not initialized properly." << std::endl;
        wcsvfree(&nwcs, &wcs_struct);  // Free if WCS failed
        return false;
    }

    std::cout << "WCS initialized successfully." << std::endl;
    // std::cout << "world at pixel 50, 50" << std::endl;
    // double pix[2] = { 50.0, 50.0 };
    // double world[2] = { 0.0, 0.0 };
    // double imgcrd[2] = { 0.0, 0.0 };
    // double phi[2] = { 0.0, 0.0 };
    // double theta[2] = { 0.0, 0.0 };
    // int stat[2] = { 0, 0 };

    // // wcsp2s now requires additional arrays (imgcrd, phi, theta) per wcslib documentation.
    // int ret = wcsp2s(wcs_ptr, 1, 2, pix, world, imgcrd, phi, theta, stat);
    // if (ret) {
    //     ERR_PRINT("wcsp2s() failed");
    // }
    // std::cout << world[0] << " " << world[1] << std::endl;

    return true;  // No need to call wcsvfree() here since we are keeping the reference
}

bool initialize_wcsprm_d(fitsfile *fptr, wcsprm **wcs_ptr) {
    int status = 0, nkeys = 0, nreject = 0, nwcs = 0;
    char *header = nullptr;
    struct wcsprm *wcs_struct = nullptr;


    if (fits_hdr2str(fptr, 1, nullptr, 0, &header, &nkeys, &status)) {
        std::cerr << "Error reading FITS header: " << status << std::endl;
        return false;
    }


    int wcspih_status = wcspih(header, nkeys, 0, 0, &nreject, &nwcs, &wcs_struct);
    free(header);  // Free header immediately after parsing
    if (wcspih_status || nwcs == 0) {
        std::cerr << "Error parsing WCS from FITS header." << std::endl;
        return false;
    }

    
    // Assign wcs_struct[0] directly to wcs_ptr, transferring ownership
    *wcs_ptr = wcs_struct;  // Now wcs_ptr points to the allocated structure

    // Finalize WCS initialization
    int wcsset_status = wcsset(*wcs_ptr);

    // Check if WCS was initialized properly
    if ((*wcs_ptr)->flag == 0) {
        std::cerr << "ERROR: WCS not initialized properly." << std::endl;
        wcsvfree(&nwcs, &wcs_struct);  // Free if WCS failed
        return false;
    }

    return true;  // No need to call wcsvfree() here since we are keeping the reference
}


bool initialize_wcsprm_old(fitsfile *fptr, wcsprm **wcs_ptr) { //struct wcsprm &wcs) {
    int status = 0, nkeys = 0, nreject = 0, nwcs = 0;
    char *header = nullptr;
    struct wcsprm *wcs_struct = nullptr;


    // Read the FITS header into a string
    if (fits_hdr2str(fptr, 1, nullptr, 0, &header, &nkeys, &status)) {
        std::cerr << "Error reading FITS header: " << status << std::endl;
        return false;
    }


    // Parse the WCS information from the header
    int wcspih_status = wcspih(header, nkeys, 0, 0, &nreject, &nwcs, &wcs_struct);
    if (wcspih_status) {
        std::cerr << "Error parsing WCS from FITS header. wcspih() returned: " << wcspih_status << std::endl;
        free(header);
        return false;
    }


    *wcs_ptr = new wcsprm(wcs_struct[0]);  // Copy first WCS
    (*wcs_ptr)->m_flag = 0;  // Prevent deallocation

    if (nwcs > 0) {
        // wcs = wcs_struct[0];
        // wcs_struct[0].m_flag = 0;  // Prevent automatic deallocation
    } else {
        std::cerr << "ERROR: No WCS structures found in the FITS header." << std::endl;
    }

    int wcsset_status = wcsset(*wcs_ptr); // &wcs);

    // Print some debug info about the WCS
    if ((*wcs_ptr)->flag) {
    } else {
        std::cerr << "ERROR: WCS not initialized properly." << std::endl;
    }

    // Free allocated memory
    // bool ok = nwcs > 0;
    wcsvfree(&nwcs, &wcs_struct);
    free(header);

    return (*wcs_ptr)->flag != 0;//ok;
}

bool FITSReader::load_fits(const String &path) {
    if (loaded && fptr) {
        fits_close_file(fptr, &status);
        loaded = false;
    }
    status = 0;
    const char *c_path = path.utf8().get_data();
    if (fits_open_file(&fptr, c_path, READONLY, &status))
        return false;
    loaded = true;
    // Move to primary HDU.
    fits_movabs_hdu(fptr, 1, nullptr, &status);

    // Retrieve image parameters (bitpix, naxis, naxes) from primary HDU.
    status = 0;
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status)) {
        naxis = 0;
    }

    status = 0;
   
   int n_hdus = 0;
fits_get_num_hdus(fptr, &n_hdus, &status);

// int hdutype = 0;
// for (int i = 1; i <= n_hdus; i++) {
//     fits_movabs_hdu(fptr, i, &hdutype, &status);
//     if (hdutype == IMAGE_HDU) {
//         std::cout << "Using HDU #" << i << " (first image HDU)" << std::endl;
//         break;
//     }
// }


    struct wcsprm wcs;
    fits_movabs_hdu(fptr, 1, nullptr, &status);

    // Initialize WCS
    if (initialize_wcsprm(fptr, &wcs_ptr)) {

    } else {
        std::cerr << "Failed to initialize WCS." << std::endl;
    }


    // wcs_ptr = &wcs;
    // Close the FITS file
    // fits_close_file(fptr, &status);

    return true;
}



PackedFloat32Array FITSReader::get_image_data() {
    PackedFloat32Array arr;
    if (!loaded || !fptr)
        return arr;
    int hdutype = 0;
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    if (fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status))
        return arr;
    const long num_pixels = naxes[0] * naxes[1];
    arr.resize(num_pixels);
    float* data_ptr = const_cast<float*>(arr.ptr());
    long fpixel = 1;
    int anynul = 0;
    status = 0;
    if (fits_read_img(fptr, TFLOAT, fpixel, num_pixels, nullptr, data_ptr, &anynul, &status))
        arr.resize(0);
    return arr;
}

Dictionary FITSReader::get_header_info() {
    Dictionary header;
    if (!loaded || !fptr)
        return header;
    int nkeys = 0;
    status = 0;
    if (fits_get_hdrspace(fptr, &nkeys, nullptr, &status))
        return header;
    for (int i = 1; i <= nkeys; i++) {
        char card[FLEN_CARD];
        if (fits_read_record(fptr, i, card, &status))
            continue;
        std::string card_str(card);
        if (card_str.find("=") == std::string::npos)
            continue;
        std::string key = card_str.substr(0, 8);
        key.erase(key.find_last_not_of(' ') + 1);
        size_t eq_pos = card_str.find("=");
        size_t val_start = eq_pos + 1;
        while (val_start < card_str.size() && card_str[val_start] == ' ')
            val_start++;
        size_t slash_pos = card_str.find("/", val_start);
        std::string value = (slash_pos != std::string::npos) ?
            card_str.substr(val_start, slash_pos - val_start) :
            card_str.substr(val_start);
        value.erase(value.find_last_not_of(' ') + 1);
        header[String(key.c_str())] = String(value.c_str());
    }
    return header;
}

Dictionary FITSReader::get_info() {
    Dictionary info;
    if (!loaded || !fptr)
        return info;
    int n_hdus = 0;
    status = 0;
    if (fits_get_num_hdus(fptr, &n_hdus, &status))
        return info;
    
    Array hdus;
    for (int i = 1; i <= n_hdus; i++) {
        Dictionary hdu_info;
        int hdutype = 0;
        status = 0;
        // Move to the absolute HDU number i:
        if (fits_movabs_hdu(fptr, i, &hdutype, &status))
            continue;
        
        hdu_info["index"] = i;
        String type_str;
        if (hdutype == IMAGE_HDU)
            type_str = "image";
        else if (hdutype == ASCII_TBL || hdutype == BINARY_TBL)
            type_str = "table";
        else
            type_str = "unknown";
        hdu_info["type"] = type_str;
        
        if (hdutype == IMAGE_HDU) {
            int img_bitpix = 0, img_naxis = 0;
            long img_naxes[10] = {0};
            if (fits_get_img_param(fptr, 10, &img_bitpix, &img_naxis, img_naxes, &status) == 0) {
                hdu_info["naxis"] = img_naxis;
                Array dims;
                for (int j = 0; j < img_naxis; j++) {
                    dims.append((int)img_naxes[j]);
                }
                hdu_info["dimensions"] = dims;
            }
        } else if (hdutype == ASCII_TBL || hdutype == BINARY_TBL) {
            long nrows = 0;
            int ncols = 0;
            if (fits_get_num_rows(fptr, &nrows, &status) == 0)
                hdu_info["nrows"] = (int)nrows;
            if (fits_get_num_cols(fptr, &ncols, &status) == 0)
                hdu_info["ncols"] = ncols;
        }
        
        char extname[FLEN_VALUE];
        status = 0;
        if (fits_read_key(fptr, TSTRING, "EXTNAME", extname, nullptr, &status) == 0)
            hdu_info["name"] = String(extname);
        
        hdus.append(hdu_info);
    }
    info["hdus"] = hdus;
    return info;
}


Dictionary FITSReader::get_table_data(int hdu_index) {
    Dictionary ret;
    if (!loaded || !fptr)
        return ret;
    status = 0;
    int hdutype = 0;
    if (fits_movabs_hdu(fptr, hdu_index, &hdutype, &status)) {
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    if (hdutype != ASCII_TBL && hdutype != BINARY_TBL) {
        printf("Table: %d\n", hdutype);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    long nrows = 0;
    int ncols = 0;
    if (fits_get_num_rows(fptr, &nrows, &status)) {
        printf("CC");
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    if (fits_get_num_cols(fptr, &ncols, &status)) {
        printf("DD");
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    Array column_names;
    Dictionary table_data;
    for (int col = 1; col <= ncols; col++) {
        char ttype[FLEN_VALUE];
        char keyname[FLEN_KEYWORD];
        std::snprintf(keyname, FLEN_KEYWORD, "TTYPE%d", col);
        status = 0;
        if (fits_read_key(fptr, TSTRING, keyname, ttype, nullptr, &status)) {
            std::snprintf(ttype, FLEN_VALUE, "COL%d", col);
            status = 0;
        }
        String colname = String(ttype);
        column_names.append(colname);
        
        std::vector<float> buffer(nrows);
        float nullval = 0.0f;
        int anynull = 0;
        if (fits_read_col(fptr, TFLOAT, col, 1, 1, nrows, &nullval, buffer.data(), &anynull, &status)) {
            ret["error"] = int(Error::FAILED);
            return ret;
        }
        PackedFloat32Array col_array;
        col_array.resize(nrows);
        for (int i = 0; i < nrows; i++) {
            col_array.set(i, buffer[i]);
        }
        table_data[colname] = col_array;
    }
    ret["error"] = int(Error::OK);
    ret["columns"] = column_names;
    ret["data"] = table_data;
    return ret;
}


Vector2 FITSReader::world_to_pixel(double ra, double dec) {
    if (!wcs_ptr) {
        ERR_PRINT("WCS not initialized");
        return Vector2();
    }
    double world[2] = { ra, dec };
    double pix[2] = { 0.0, 0.0 };
    double phi[2] = { 0.0, 0.0 };
    double theta[2] = { 0.0, 0.0 };
    double imgcrd[2] = { 0.0, 0.0 };
    int stat[2] = { 0, 0 };

    // wcss2p now requires additional arrays (phi, theta, imgcrd) per wcslib documentation.
    int ret = wcss2p(wcs_ptr, 1, 2, world, phi, theta, imgcrd, pix, stat);

    if (ret) {
        ERR_PRINT("wcss2p() failed");
        return Vector2();
    }
    return Vector2(pix[0], pix[1]);
}

Vector2 FITSReader::pixel_to_world(double x, double y) {
    if (!wcs_ptr) {
        ERR_PRINT("WCS not initialized");
        return Vector2();
    }

    double pix[2] = { x, y };
    double world[2] = { 0.0, 0.0 };
    double imgcrd[2] = { 0.0, 0.0 };
    double phi[2] = { 0.0, 0.0 };
    double theta[2] = { 0.0, 0.0 };
    int stat[2] = { 0, 0 };


    int ret = wcsp2s(wcs_ptr, 1, 2, pix, imgcrd, phi, theta, world, stat);

    
    if (ret) {
        ERR_PRINT("wcsp2s() failed");
        return Vector2();
    }
    
    return Vector2(world[0], world[1]);
}


// Vector2 FITSReader::pixel_to_world_old(double x, double y) {
//     if (!wcs_ptr) {
//         ERR_PRINT("WCS not initialized");
//         return Vector2();
//     }
//     double pix[2] = { x, y };
//     double world[2] = { 0.0, 0.0 };
//     double imgcrd[2] = { 0.0, 0.0 };
//     double phi[2] = { 0.0, 0.0 };
//     double theta[2] = { 0.0, 0.0 };
//     int stat[2] = { 0, 0 };

//     // wcsp2s now requires additional arrays (imgcrd, phi, theta) per wcslib documentation.
//     int ret = wcsp2s(wcs_ptr, 1, 2, pix, world, imgcrd, phi, theta, stat);
//     if (ret) {
//         ERR_PRINT("wcsp2s() failed");
//         return Vector2();
//     }
//     return Vector2(world[0], world[1]);
// }


Dictionary FITSReader::read_spectrum(const String &filename, const String &column_name) {
    Dictionary ret;
    fitsfile *fptr_local = nullptr;
    int status_local = 0;
    long nrows = 0;
    int colnum = 0;

    if (fits_open_file(&fptr_local, filename.utf8().get_data(), READONLY, &status_local)) {
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    int hdutype;
    if (fits_movabs_hdu(fptr_local, 2, &hdutype, &status_local)) {
        fits_close_file(fptr_local, &status_local);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    if (hdutype != ASCII_TBL && hdutype != BINARY_TBL) {
        fits_close_file(fptr_local, &status_local);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    if (fits_get_num_rows(fptr_local, &nrows, &status_local)) {
        fits_close_file(fptr_local, &status_local);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    std::string colname_std = column_name.utf8().get_data();
    if (fits_get_colnum(fptr_local, CASEINSEN, &colname_std[0], &colnum, &status_local)) {
        fits_close_file(fptr_local, &status_local);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    std::vector<float> buffer(nrows);
    float nullval = 0.0f;
    int anynull = 0;
    if (fits_read_col(fptr_local, TFLOAT, colnum, 1, 1, nrows, &nullval, buffer.data(), &anynull, &status_local)) {
        fits_close_file(fptr_local, &status_local);
        ret["error"] = int(Error::FAILED);
        return ret;
    }
    fits_close_file(fptr_local, &status_local);
    PackedFloat32Array spectrum_data;
    spectrum_data.resize(nrows);
    for (int i = 0; i < nrows; i++) {
        spectrum_data.set(i, buffer[i]);
    }
    ret["error"] = int(Error::OK);
    ret["spectrum_data"] = spectrum_data;
    return ret;
}

