/**
 * Copyright © 2016 Andrey Nevolin, https://github.com/AndreyNevolin
 * Twitter: @Andrey_Nevolin
 * LinkedIn: https://www.linkedin.com/in/andrey-nevolin-76387328
 *
 * Reference implementation of "split_FindBound()" function intended
 * for splitting of FASTA files (FASTA is a bioinformatics format for
 * representing nucleotide or peptide sequences). The function
 * recognizes the simplest version of the format:
 * ">IDENTIFIER\n"
 * "SEQUENCE\n"
 *
 * To split files of a different format, one needs to replace the
 * contents of the function with appropriate code
 */

/**
 * Short description of supported file format. This value will be
 * shown in a help message explaining the purpose of the tool
 */
#define SPLIT_FILE_FORMAT_NAME "FASTA"

/**
 * Find actual bound of an element inside a buffer. The actual bound
 * should preferably be close to a projected bound
 *
 * Input: buff - pointer to the buffer
 *        projected_bound - desired offset of an element's last byte
 *        buff_size - buffer size
 *        is_first_block - indicator that a new output file will be started
 *                         with the data in the buffer (i.e. we're going
 *                         to write first chunk of the file)
 *
 * Return value: offset relative to "buff" if element bound was found, or
 *               SPLIT_BOUND_NOT_FOUND in the opposite case.
 *               "offset" should be from "-1" to "buff_size - 1".
 *               "-1" means that first byte of "buff" concides with first
 *               byte of some element and we want to finish current output
 *               file without adding any "buff" contents to it. This implies
 *               that the output file already contains an integer number of
 *               "elements". This is only possible when we don't start a
 *               new output file but continue adding data to already
 *               existing one. So, returning "-1" is only allowed when
 *               "is_first_block" indicator is "false". This is exactly the
 *               reason of having this indicator among the parameters: to
 *               allow to not add more data to an output file. Consider an
 *               example: buffer contents starts with a new element of
 *               size 100, projected bound is 10. The choice is adding 89
 *               to the desired size vs. subtracting 11. Subtracting 11
 *               might result in a better balancing of output file sizes
 */
int64_t split_FindBound( const char * const buff,
                         int64_t projected_bound,
                         int64_t buff_size,
                         bool is_first_block)
{
    /* Check that the desired bound falls inside the buffer */
    SPLIT_ASSERT( projected_bound >= 0);
    SPLIT_ASSERT( projected_bound < buff_size);

    int64_t distance_to_l_bound = projected_bound + 1;
    int64_t distance_to_u_bound = buff_size - projected_bound;
    char num_new_lines = 0, right_symb = 0, left_symb = 0;

    /* Seek from the desired bound byte-by-byte to left and right simultaneously.
       Stop when we find an element bound. This bound is guaranteed to be the
       closest to the desired bound. We check left byte before the right one,
       because it's preferrable for us to have output files that are smaller than
       projected. That's becasue the last output file is generally expected to be
       smaller than other pieces. Keeping other pieces smaller might result in a
       bigger last file and - hence - in better size balancing (though it's still
       just a heuristic) */
    for ( int64_t i = 0; i < std::max( distance_to_l_bound, distance_to_u_bound); i++ )
    {
        if ( i >= distance_to_l_bound )
        {
            /* Shouldn't seek left to the left bound. So, use arbitrary symbol
               which is not equal to the element start symbol */
            left_symb = 'a';
        } else
        {
            left_symb = buff[projected_bound - i];
        }

        if ( i >= distance_to_u_bound )
        {
            /* Shouldn't seek right to the right bound. So, use arbitrary symbol
               which is not equal to the element start symbol */
            right_symb = 'a';
        } else
        {
            right_symb = buff[projected_bound + i];
        }

        /* Element start symbol found */
        if ( left_symb == '>' )
        {
            if ( (i < distance_to_l_bound - 1) || !is_first_block )
            {
                /* Return previous symbol only in case when current buffer
                   contents is not intended to start new output file */
                return projected_bound - i - 1;
            }
        }

        /* Element start symbol found */
        if ( (right_symb == '>') && (i != 0) )
        {
            /* Return previous symbol only in case when it doesn't coinside
               with the projected bound (if it does that means that this
               bound doesn't work for us, because otherwise we wouldn't
               reach this condition. The condition on "left_symb" above
               would "return" in this case) */
            return projected_bound + i - 1;
        }

        if ( left_symb == '\n' )
        {
            num_new_lines++;
        }

        if ( (right_symb == '\n') && (i != 0) )
        {
            num_new_lines++;
        }

        /* Each element should have exactly two newlines. So, if we detected
           two newlines, the 'upper' one may be used as an element bound. But
           here we use a stronger condition: we reached upper bound of the
           buffer and it represents the second newline (i.e. last byte of the
           buffer conicides with an end of some element */
        if ( (num_new_lines == 2) && (i >= distance_to_u_bound - 1)
             && (buff[buff_size - 1] == '\n') )
        {
            return buff_size - 1;
        }

        SPLIT_ASSERT( num_new_lines <= 2);
    }

    return SPLIT_BOUND_NOT_FOUND;
}
