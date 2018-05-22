/**
 * Copyright © 2016 Andrey Nevolin, https://github.com/AndreyNevolin
 * Twitter: @Andrey_Nevolin
 * LinkedIn: https://www.linkedin.com/in/andrey-nevolin-76387328
 *
 * Utility that splits a file of "elements" into pieces of roughly equivalent sizes.
 * Each piece would contain an integer number of "elements" (i.e. one single element
 * couldn't be divided between pieces. Each element is required to belong to some
 * piece as a whole).
 *
 * Exact definition of "elements" (and thus complete purpose of the tool) is provided
 * by a party that builds the tool.
 *
 * The tool is incomplete. It's a prototype which is agnostic of exact source file
 * format (i.e. agnostic of "elements" definition). To make the tool complete, one
 * needs to implement "split_FindBound()" function which is defined in
 * "find_bound.cpp" file. This file comes with a reference implementation of
 * "split_FindBound()" intended for splitting files in FASTA format (bioinformatics
 * format for representing nucleotide or peptide sequences). To build the tool for
 * splitting files of different format, one needs to replace the contents of
 * "split_FindBound()" with corresponding code.
 *
 * The tool reads data from input file and writes data to output files in aligned
 * chunks of fixed size (except maybe the last chunk read from/appended to a file).
 * The size of a chunk can be specified by a user. By default it's 4Mb.
 *
 * To ensure that output files contain an integer number of "elements" each, the tool
 * somehow needs to recognize bounds of individual elements. When last chunk of data
 * is added to an output file, projected bound of a file might be shifted up or down
 * to make the file include integer number of elements. That's where "split_FindBound()"
 * comes to play. Given a buffer with data and projected output file bound, it should
 * be able to find an element bound which is close to the projected bound.
 *
 * The tool works best when the chunk size is much bigger than the size of any
 * "element". For more information about "split_FindBound()" see comments inside
 * "find_bound.cpp"
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#ifdef SPLIT_DEBUG
#include <execinfo.h>
#endif
#include <string>

/* Default buffer size is 4Mb */
#define SPLIT_BUFFER_SIZE_DEFAULT 4194304

#ifdef SPLIT_DEBUG
/**
 * Print stack dump
 */
static inline void split_PrintStack()
{
    void *frames[100];
    size_t stack_depth;
    char **calls; 

    fprintf( stderr, "Stack trace: \n");
    stack_depth = backtrace( frames, 100);
    calls = backtrace_symbols( frames, stack_depth);

    for ( size_t i = 0; i < stack_depth; i++ )
    {
        printf ("\t[%2lu] %s\n", stack_depth - i - 1, calls[i]);
    }

    free( calls); 
}

#define SPLIT_ABORT \
            fprintf( stderr, "\n"), \
            fprintf( stderr, "Internal error: file \"%s\", line %u\n", __FILE__, __LINE__), \
            fprintf( stderr, "\n"), \
            fflush( NULL), \
            split_PrintStack(), \
            abort()
#endif /* SPLIT_DEBUG */

#ifdef SPLIT_DEBUG
#    define SPLIT_ASSERT( condition_) ((condition_) || ( SPLIT_ABORT, 0))
#else
#    define SPLIT_ASSERT( condition_)
#endif

#define SPLIT_OUT( format_, ...) \
            fprintf( stdout, format_, ##__VA_ARGS__); \
            fprintf( stdout, "\n");

#define SPLIT_ERROR( format_, ...) \
            fprintf( stderr, format_, ##__VA_ARGS__); \
            fprintf( stderr, "\n"); \
            exit( EXIT_FAILURE);

/**
 * Wrapper around "strerror_r()" aimed at making calls to
 * "strerror_r()" more compact and manageable
 */
static inline char *SPLIT_STRERROR_R( char *buff, size_t size)
{
    SPLIT_ASSERT( buff);

    return strerror_r( errno, buff, size);
}

/* Value that should be returned by "split_FoundBound()" function if
   element bound wasn't found */
#define SPLIT_BOUND_NOT_FOUND -12345
#include "find_bound.cpp"

/**
 * Structure to keep command-line and derived options
 */
typedef struct
{
    /* Path to input file */
    std::string input_path;
    /* Name of output directory */
    std::string output_dir;
    /* Base name of output files */
    std::string output_file;
    /* Number of output files */
    int64_t num_pieces;
    /* Size in bytes of a buffer used to read/write files. Data
       will be read/written mostly in chunks of this size */
    int64_t buffer_size;
} split_Opts_t;

/**
 * Description of command line options used by "getopt_long()"
 */
static struct option split_long_options_desc[] =
{
    /* Output directory */
    {"od", required_argument, 0, 'd'},
    /* Base of output file name */
    {"of", required_argument, 0, 'f'},
    /* Chunk size */
    {"cs", required_argument, 0, 'c'},
    {0,    0,                 0, 0}
};

/**
 * Description of short options
 * 
 * '+' means that parsing stops when first non-option argument is encountered
 * ':' in the beginning means that ':' will be returned by 'getopt()' for missing
 *     option arguments
 * ':' after option character means that this particular option requires an
 *     argument
 * 'n' - number of output files
 */
std::string split_short_options_desc = "+:n:";

/**
 * Print simple assistance message and exit
 */
void split_ExitWithAssist( const char *buff, const char *prog_name)
{
    SPLIT_ERROR( "%s\n\nRun \"%s\" without options for usage guidelines",
                 buff, prog_name);
}

static const char *usage_format[] =
{
    "Usage: %s -n <number of pieces> [-od <output directory>] "
    "[-of <basis for output file name>] [-cs <chunk size>] <path to file to split>",
    ""
};

static const char *usage_options[] =
{
    " ",
    "Split file in " SPLIT_FILE_FORMAT_NAME " format into pieces of roughly equivalent "
    "sizes. Each piece will contain an integer number of records from the input file",
    " ",
    "OPTIONS:",
    "   -n          Number of pieces to produce. Each piece will be placed into",
    "               a separate file named \"<file name>.<number>\", where \"<file",
    "               name>\" is a name provided through \"--of\" option or the",
    "               name of the input file (if \"--of\" option is not used).",
    "               \"<number>\" is a sequential number of a piece",
    "       --od    Path to output directory. By default current directory will",
    "               be used for output",
    "       --of    Basis for output file names. Output files will be named",
    "               \"<file name>.<number>\", where \"<file name>\" is a string",
    "               provided through this option or the name of the input file",
    "               (if this option is not used)",
    "       --cs    Data chunk size. Data will be read from input file and written",
    "               to output files by chunks of this size (except maybe the last",
    "               chunk read from/written to file). The size may be provided in",
    "               different units: 512B, 4K, 8M, 1G (\"B\" for bytes, \"K\" for",
    "               kilobytes, \"M\" for megabytes, and \"G\" for gigabytes). If",
    "               units identifier is omitted, byte units are implied. The",
    "               default value for this option is 4M",
    ""
};

/**
 * Print help message for this program
 */
void split_PrintHelp( const char *const prog_name)
{
    for ( int i = 0; ; i++ )
    {
        if ( !strlen( usage_format[i]) )
        {
            break;
        }

        SPLIT_OUT( usage_format[i], prog_name);
    }

    for ( int i = 0; ; i++ )
    {
        if ( !strlen( usage_options[i]) )
        {
            break;
        }

        SPLIT_OUT( "%s", usage_options[i]);
    }
}

/**
 * Check if 'strtol' conversion was successful
 */
bool split_IsStrtolOK( char first_symb, int err_no, char curr_symb, char base)
{
    SPLIT_ASSERT( (base > 1) && (base < 11));

    /* Here we don't allow leading '+' or '-' which is normally
       allowed by 'strtol' */
    bool check_res = (first_symb < '0')
                     || (first_symb > '0' + base - 1)
                     || (err_no == ERANGE)
                     || (curr_symb != '\0');

    return !check_res;
}

/**
 * Initialize options
 */
int split_InitOpts( split_Opts_t *opts)
{
    opts->buffer_size = SPLIT_BUFFER_SIZE_DEFAULT;
    opts->num_pieces = 0;

    return 0;
}

/**
 * Parse command line
 */
int split_ParseCmdLine( int argc, char *argv[], split_Opts_t *opts)
{
    int long_opt_index = -1;
    int getopt_res = -1;
    std::string prog_name = basename( argv[0]);

    while ( 1 )
    {
        int curr_optind = optind;

        getopt_res = getopt_long( argc, argv, split_short_options_desc.c_str(),
                                  split_long_options_desc, &long_opt_index);

        /* This return code means the input was fully parsed */
        if ( getopt_res == -1 )
        {
            break;
        }

        /* Buffer for error messages */
        char buff[200];

        switch ( getopt_res )
        {
            /* Output directory */
            case 'd':
                opts->output_dir = std::string( optarg);

                break;

            /* Output file base name */
            case 'f':
                opts->output_file = std::string( optarg);

                break;

            /* Chunk size */
            case 'c':
            {
                std::string arg_copy( optarg);
                char unit = 0;

                /* Check if units indentifier was provided */
                if ( arg_copy.back() < '0' || arg_copy.back() > '9' )
                {
                    /* Save units identifier */
                    unit = arg_copy.back();
                    /* Delete units identifier from the string */
                    arg_copy.erase( arg_copy.size() - 1, 1);
                }

                char *c_ptr = 0;
                int64_t buff_size = strtol( arg_copy.c_str(), &c_ptr, 10);

                if ( !split_IsStrtolOK( arg_copy.front(), errno, *c_ptr, 10) )
                {
                    split_ExitWithAssist( "Integer with units is expected "
                                          "for chunk size", prog_name.c_str());
                }

                int shift_val = 0;

                switch ( unit )
                {
                    /* Units is byte or units wasn't provided */
                    case 0:
                    case 'b':
                    case 'B':
                        break;

                    /* Kilobytes */
                    case 'k':
                    case 'K':
                        shift_val = 10;

                        break;

                    /* Megabytes */
                    case 'm':
                    case 'M':
                        shift_val = 20;

                        break;

                    /* Gigabytes */
                    case 'g':
                    case 'G':
                        shift_val = 30;

                        break;

                    default:
                        split_ExitWithAssist( "Unexpected units identifier for buffer "
                                              "size", prog_name.c_str());
                }

                /* We add "1" to "shift_val" below because the programm uses internally
                   a buffer of double size */
                if ( buff_size > (INT64_MAX >> (shift_val + 1)) )
                {
                    SPLIT_ERROR( "Chunk size if too big. Maximum size is %ld bytes",
                                 INT64_MAX / 2);
                }

                opts->buffer_size = buff_size << shift_val;
 
                break;
            }

            /* Number of pieces */
            case 'n':
            {
                char *c_ptr = 0;

                opts->num_pieces = strtol( optarg, &c_ptr, 10);

                if ( !split_IsStrtolOK( optarg[0], errno, *c_ptr, 10) )
                {
                    split_ExitWithAssist( "Integer is expected for number of pieces",
                                          prog_name.c_str());
                }

                if ( opts->num_pieces < 2 )
                {
                    SPLIT_ERROR( "Number of pieces should be greater than 1");
                }

                break;
            }

            /* Missing mandatory argument */
            case ':':
                snprintf( buff, sizeof( buff),
                         "Mandatory argument is missing for \"%s\"",
                         argv[curr_optind]);
                split_ExitWithAssist( buff, prog_name.c_str());

                break;

            /* Unknow option */
            case '?':
                snprintf( buff, sizeof( buff), "Unknown option: %s", argv[curr_optind]);
                split_ExitWithAssist( buff, prog_name.c_str());

                break;

            /* Unexpected return value of 'getopt' */
            default:
                SPLIT_ASSERT( 0);
        }
    }

    if ( optind == argc )
    {
        if ( argc == 1 )
        {
            split_PrintHelp( prog_name.c_str());

            exit( EXIT_SUCCESS);
        } else
        {
            split_ExitWithAssist( "Name of input file is required", prog_name.c_str());
        }
    }

    if ( !opts->num_pieces )
    {
        split_ExitWithAssist( "Number of pieces is required", prog_name.c_str());
    }

    opts->input_path = std::string( argv[optind]);

    if ( optind < argc - 1 )
    {
        SPLIT_OUT( "Warning: several input file names were provided. Only first one "
                   "will be used");
    }

    if ( (opts->output_file).empty() )
    {
        char *buff = (char *)malloc( opts->input_path.size() + 1);

        if ( !buff )
        {
            SPLIT_ERROR( "Couldn't allocate memory for a temporary buffer");
        }

        strcpy( buff, opts->input_path.c_str());
        /* Use input file name as the base name for pieces */
        opts->output_file = basename( buff);
        free( buff);
    }

    if ( (opts->output_dir).empty() )
    {
        /* Use current directory for output */
        opts->output_dir = ".";
    }

    return 0;
}

/**
 * Calculate number of decimal digits needed to write down N-1
 */
static int split_CalcNumWidth( int64_t num)
{
    int num_digits = 0;

    SPLIT_ASSERT( num > 1);
    num--;

    while ( num )
    {
        num_digits++;
        num /= 10;
    }

    return num_digits;
}

/**
 * Get size of an input file
 */
int64_t split_GetInputSize( int fd)
{
    /* Check that the file is at zero offset. Otherwise we would need
       to save current offset to return to it later */
    SPLIT_ASSERT( !lseek( fd, 0, SEEK_CUR));

    char err_msg[500];
    int64_t input_size = lseek( fd, 0, SEEK_END);

    if ( input_size == -1 )
    {
        SPLIT_ERROR( "Cannot seek input file: %s",
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    }

    /* Return to the beginning of the file */
    if ( lseek( fd, 0, SEEK_SET) == -1 )
    {
        SPLIT_ERROR( "Cannot seek input file: %s",
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    }

    return input_size;
}

/**
 * Create output file for a new piece
 */
static int split_StartNewPiece( const std::string & dir_name,
                                const std::string & file_name,
                                int num_digits,
                                int64_t piece_num)
{
    std::string output_name = dir_name + "/" + file_name + ".";
    char buff[64];
    char err_msg[500];

    snprintf( buff, sizeof( buff), "%0*ld", num_digits, piece_num);
    output_name += buff;

    int output_fd = open( output_name.c_str(), O_CREAT | O_EXCL | O_WRONLY);

    if ( output_fd == -1 )
    {
        SPLIT_ERROR( "Cannot create output file \"%s\": %s", output_name.c_str(),
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    }

    return output_fd;
}

/**
 * Sync output file to persistent store and close
 */
static int split_FinalizePiece( int fd, int64_t piece_num)
{
    char err_msg[500];
    int64_t piece_size = lseek( fd, 0, SEEK_END);

    if ( fsync( fd) == -1 )
    {
        SPLIT_ERROR( "Cannot sync output file: %s",
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    }

    close( fd);

    char prefix[100];

    snprintf( prefix, sizeof( prefix), "Piece %ld written. Size: ", piece_num + 1);

    if ( piece_size != -1 )
    {
        char units = 0;
        double value = 0;

        if ( piece_size / (1024 * 1024 * 1024) )
        {
            value = piece_size / (1024 * 1024 * 1024);
            units = 'G';
        } else if ( piece_size / (1024 * 1024) )
        {
            value = piece_size / (1024 * 1024);
            units = 'M';
        } else if ( piece_size / 1024 )
        {
            value = piece_size / 1024;
            units = 'K';
        }

        if ( units )
        {
            SPLIT_OUT( "%s%.1f%c (%ld bytes)", prefix, value, units, piece_size);
        } else
        {
            SPLIT_OUT( "%s%ld bytes", prefix, piece_size);
        }
    } else
    {
        SPLIT_OUT( "%sunknown", prefix);
    }

    return 0;
}

/**
 * Read chunk of data from an input file to upper half of the double-buffer
 */
static int64_t split_FillUpperBuffHalfFromInput( int fd,
                                                 char *double_buff,
                                                 int64_t buff_size,
                                                 int64_t bytes_available)
{
    if ( !bytes_available )
    {
        /* End of input file was reached */
        return 0;
    }

    int64_t io_size = buff_size;

    /* We always read data in chunks of fixed size, except maybe the last
       read of input file (because size of input file may not be multiple
       of chunk size). So, even if amount of data needed to complete current
       piece is less than chunk size, we still read full chunk. If we don't
       add it fully to the current piece, remaining data will be used to
       start next piece */
    if ( io_size > bytes_available )
    {
        io_size = bytes_available;
    }

    char err_msg[500];
    int64_t bytes_read = read( fd, double_buff + buff_size, io_size);

    if ( bytes_read == -1 )
    {
        SPLIT_ERROR( "Cannot read data from the input file: %s",
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    } else if ( bytes_read != io_size )
    {
        SPLIT_ERROR( "Read %ld bytes from the input file. %ld bytes were expected. "
                     "Is it a regular file?", bytes_read, io_size);
    }

    return bytes_read;
}

/**
 * Write chunk of data to output file
 */
int64_t split_WriteOutput( int fd, char *buff, int64_t data_start, int64_t data_end)
{
    char err_msg[500];
    int64_t io_size = data_end - data_start + 1;
    int64_t bytes_written = write( fd, buff + data_start, io_size);

    if ( bytes_written == -1 )
    {
        SPLIT_ERROR( "Cannot write data to output file: %s",
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    } else if ( bytes_written != io_size )
    {
        SPLIT_ERROR( "Written %ld bytes to an output file. %ld bytes were expected. "
                     "Is it a regular storage device?", bytes_written, io_size);
    }

    return bytes_written;
}

/**
 * Determine upper bound of data-chunk that will be transferred from
 * the double-buffer to output file
 */
int64_t split_CalcUpperBoundOfOutputTransfer( char *buff,
                                              int64_t buff_size,
                                              int64_t data_start,
                                              int64_t data_end,
                                              int64_t projected_max,
                                              bool is_first_block,
                                              bool is_end_of_input,
                                              bool is_last_piece)
{
    /* Size of active data currently in the double-buffer */
    int64_t active_data_size = data_end - data_start + 1;
    int64_t bound = -1;

    SPLIT_ASSERT( active_data_size < 2 * buff_size);

    /* Check if we need more data than there is available in the buffer */
    if ( projected_max > active_data_size )
    {
        /* Check if we have more than a full chunk */
        if ( active_data_size >= buff_size )
        {
            return data_start + buff_size - 1;
        } else
        {
            return data_end;
        }
    } else
    {
        /* Do not search for a closest element bound if we're writing the last
           piece. We have no choice but to append all remaining buffer contents
           to this file */
        if ( is_last_piece )
        {
            SPLIT_ASSERT( projected_max == active_data_size);

            return data_end;
        }

        /* Find element bound which is closest to projected file end */
        bound = split_FindBound( buff + data_start, projected_max - 1,
                                 data_end - data_start + 1, is_first_block);

        if ( bound != SPLIT_BOUND_NOT_FOUND )
        {
            SPLIT_ASSERT( (bound >= -1) && (bound < data_end - data_start + 1));
            SPLIT_ASSERT( (bound != -1) || !is_first_block);

            return bound + data_start;
        } else if ( is_end_of_input )
        {
            /* No more data remains in input file. Best we can do
               is to append all active data to the current piece */
            return data_end;
        } else
        {
            SPLIT_ASSERT( active_data_size >= buff_size);
            SPLIT_ERROR( "No item bound found inside a data chunk. Buffer size should be "
                         "bigger than size of any item");
        }
    }

    /* Shouldn't be here */
    SPLIT_ASSERT( 0);

    /* Return statement just to make compiler happy */
    return -1;
}

/**
 * Split source file into pieces
 *
 * This is a routine that does all significant work
 *
 * The algorithm reads and writes data in chunks of fixed size (except maybe
 * last read of input file and last write of output file). The algorithm uses
 * a buffer of double chunk size. The following rules apply to this buffer:
 * 1) when data is being read from input file it first comes to this buffer.
 *    Data contained in the buffer is called 'active data'
 * 2) first chunk of data fetched from input file comes to the second half of
 *    the buffer
 * 3) then some piece of acitve data (starting from left bound of these data)
 *    is moved to an output file
 * 4) if after that some active data remains in the buffer, it is moved to
 *    the first half of the buffer, so that right bound of the data conicides
 *    with the right bound of the first half
 * 5) the second half is filled from the input file. After that active data still
 *    represents continuous piece of input data (because when moving data from
 *    the second half to the first, we joined it to the beginning of the second
 *    half)
 * 6) steps 3)-5) are repeated, with addition that in step 4) the data is moved
 *    from the second half to the first only when left bound of active data
 *    becomes bigger (in terms of offset) than left bound of the second half
 */
int split_SplitSource( const split_Opts_t* const opts)
{
    char err_msg[500];
    int fd_input = -1;
    int64_t buff_size = opts->buffer_size;

    /* Open input file */
    fd_input = open( (opts->input_path).c_str(), O_RDONLY);

    if ( fd_input == -1 )
    {
        SPLIT_ERROR( "Cannot open file \"%s\": %s", (opts->input_path).c_str(),
                     SPLIT_STRERROR_R( err_msg, sizeof( err_msg)));
    }

    /* Calculate number of digits needed to write down number of pieces.
       We need this number because we are going to append ordinal numbers of
       pieces to their names */
    int num_digits = split_CalcNumWidth( opts->num_pieces);
    /* Get size of input file */
    int64_t input_size = split_GetInputSize( fd_input);
    int64_t bytes_available = input_size, bytes_not_read = input_size;

    /* Allocate double-buffer */
    SPLIT_ASSERT( buff_size <= (INT64_MAX / 2));

    char *double_buff = (char *)malloc( 2 * buff_size);

    if ( !double_buff )
    {
        SPLIT_ERROR( "Couldn't allocate internal buffer of size %ld",
                     2 * buff_size);
    }

    /* Initialize bounds of active data */
    int64_t data_start = buff_size;
    int64_t data_end = data_start - 1;

    for ( int64_t piece_num = 0; piece_num < opts->num_pieces; piece_num++ )
    {
        /* Calculate projected size of current piece */
        /* Divide remaining data equally between remaining pieces */
        int64_t to_read = bytes_available / (opts->num_pieces - piece_num);

        if ( bytes_available % (opts->num_pieces - piece_num) )
        {
            to_read++;
        }

        if ( !to_read )
        {
            SPLIT_ERROR( "Couldn't produce the requested number of pieces. "
                         "Only %ld pieces were writted", piece_num);
        }

        /* Start new piece */
        int output_fd = split_StartNewPiece( opts->output_dir, opts->output_file,
                                             num_digits, piece_num);
        int is_first_block = true;

        while ( to_read )
        {
            SPLIT_ASSERT( data_end >= buff_size - 1);
            SPLIT_ASSERT( data_start <= buff_size);

            int64_t bytes_read = 0;

            if ( data_end == buff_size - 1 )
            {
                bytes_read = split_FillUpperBuffHalfFromInput( fd_input, double_buff,
                                                               buff_size,
                                                               bytes_not_read);
                data_end += bytes_read;
                bytes_not_read -= bytes_read;
            }

            /* Calculate upper bound of data that will be written to output file */
            int64_t output_chunk_end = -1;
            bool is_last_piece = (piece_num == opts->num_pieces - 1);

            output_chunk_end = split_CalcUpperBoundOfOutputTransfer( double_buff,
                                                                     buff_size,
                                                                     data_start,
                                                                     data_end,
                                                                     to_read,
                                                                     is_first_block,
                                                                     !bytes_not_read,
                                                                     is_last_piece);

            if ( output_chunk_end < data_start )
            {
                SPLIT_ASSERT( output_chunk_end == data_start - 1);
                to_read = 0;
            } else if ( output_chunk_end - data_start + 1 > to_read )
            {
                to_read = 0;
            } else
            {
                to_read -= output_chunk_end - data_start + 1;
            }

            /* Append the chunk to the current output piece */
            split_WriteOutput( output_fd, double_buff, data_start, output_chunk_end);
            bytes_available -= output_chunk_end - data_start + 1;
            /* Shift left bound of active data */
            SPLIT_ASSERT( output_chunk_end < INT64_MAX);
            data_start = output_chunk_end + 1;
            is_first_block = false;

            /* Move data from upper half of the double-buffer to lower half */
            if ( (data_start >= buff_size) && (data_end - data_start >= 0) )
            {
                int64_t active_data_size = data_end - data_start + 1;

                memcpy( double_buff + buff_size - active_data_size,
                        double_buff + data_start, active_data_size);
                data_start = buff_size - active_data_size;
                data_end = buff_size - 1;
            }

            /* No data remains in the double-buffer. Need to reinitialize bounds
               of active data */
            if ( data_start > data_end )
            {
                SPLIT_ASSERT( (data_start == buff_size * 2)
                              || (!bytes_available && !bytes_not_read));
                SPLIT_ASSERT( data_end == data_start - 1);
                data_start = buff_size;
                data_end = data_start - 1;
            }
        }

        split_FinalizePiece( output_fd, piece_num);
    }

    SPLIT_ASSERT( bytes_available == 0);

    return 0;
}

int main( int argc, char *argv[])
{
    split_Opts_t opts;

    split_InitOpts( &opts);
    split_ParseCmdLine( argc, argv, &opts);
    split_SplitSource( &opts);

    exit( EXIT_SUCCESS);
}
