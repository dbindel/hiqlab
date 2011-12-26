/* Autogenerated by mclmcr/genAPI.pl: Tue Jul 26 23:11:49 2005. */

#ifndef mclcppAPI_h
#define mclcppAPI_h

/* Forward declarations */
class ref_count_obj;
class char_buffer;
class array_ref;
class array_buffer;
class error_info;

/* Class declarations */

class ref_count_obj
{
public:
    virtual int addref() = 0;
    virtual int release() = 0;
};

class char_buffer: public ref_count_obj
{
public:
    virtual int size() = 0;
    virtual const char* get_buffer() = 0;
    virtual int set_buffer(const char* str) = 0;
    virtual int compare_to(char_buffer* p) = 0;
};

class array_ref: public ref_count_obj
{
public:
    virtual mxClassID classID() = 0;
    virtual array_ref* deep_copy() = 0;
    virtual array_ref* shared_copy() = 0;
    virtual array_ref* serialize() = 0;
    virtual int element_size() = 0;
    virtual int number_of_elements() = 0;
    virtual int number_of_nonzeros() = 0;
    virtual int maximum_nonzeros() = 0;
    virtual int number_of_dimensions() = 0;
    virtual array_ref* get_dimensions() = 0;
    virtual int number_of_fields() = 0;
    virtual char_buffer* get_field_name(int i) = 0;
    virtual bool is_empty() = 0;
    virtual bool is_sparse() = 0;
    virtual bool is_numeric() = 0;
    virtual bool is_complex() = 0;
    virtual int make_complex() = 0;
    virtual bool equals(array_ref* p) = 0;
    virtual int compare_to(array_ref* p) = 0;
    virtual int hash_code() = 0;
    virtual char_buffer* to_string() = 0;
    virtual array_ref* row_index() = 0;
    virtual array_ref* column_index() = 0;
    virtual array_ref* get(int num_indices, const int* index) = 0;
    virtual array_ref* get(const char* name, int num_indices, const int* index) = 0;
    virtual array_ref* getV(int num_indices, va_list vargs) = 0;
    virtual array_ref* getV(const char* name, int num_indices, va_list vargs) = 0;
    virtual int set(array_ref* p) = 0;
    virtual array_ref* real() = 0;
    virtual array_ref* imag() = 0;
    virtual int get_numeric(mxDouble* x, int len) = 0;
    virtual int get_numeric(mxSingle* x, int len) = 0;
    virtual int get_numeric(mxInt8* x, int len) = 0;
    virtual int get_numeric(mxUint8* x, int len) = 0;
    virtual int get_numeric(mxInt16* x, int len) = 0;
    virtual int get_numeric(mxUint16* x, int len) = 0;
    virtual int get_numeric(mxInt32* x, int len) = 0;
    virtual int get_numeric(mxUint32* x, int len) = 0;
    virtual int get_numeric(mxInt64* x, int len) = 0;
    virtual int get_numeric(mxUint64* x, int len) = 0;
    virtual int get_char(mxChar* x, int len) = 0;
    virtual int get_logical(mxLogical* x, int len) = 0;
    virtual int set_numeric(const mxDouble* x, int len) = 0;
    virtual int set_numeric(const mxSingle* x, int len) = 0;
    virtual int set_numeric(const mxInt8* x, int len) = 0;
    virtual int set_numeric(const mxUint8* x, int len) = 0;
    virtual int set_numeric(const mxInt16* x, int len) = 0;
    virtual int set_numeric(const mxUint16* x, int len) = 0;
    virtual int set_numeric(const mxInt32* x, int len) = 0;
    virtual int set_numeric(const mxUint32* x, int len) = 0;
    virtual int set_numeric(const mxInt64* x, int len) = 0;
    virtual int set_numeric(const mxUint64* x, int len) = 0;
    virtual int set_char(const mxChar* x, int len) = 0;
    virtual int set_logical(const mxLogical* x, int len) = 0;
};

class array_buffer: public ref_count_obj
{
public:
    virtual int size() = 0;
    virtual array_ref* get(int offset) = 0;
    virtual int set(int offset, array_ref* p) = 0;
    virtual int add(array_ref* pa) = 0;
    virtual int remove(int offset) = 0;
    virtual int clear() = 0;
    virtual array_ref* to_cell(int offset, int len) = 0;
};

class error_info: public ref_count_obj
{
public:
    virtual const char* get_message() = 0;
};

/* Multiple include guard */
#endif
