/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*                                                  												*/
/*         							  Copyright © 2019 Dominik Řičánek								*/
/*											All rights reserved             						*/
/****************************************************************************************************/

/**
 * @file exception_handler.h
 * @brief Exception handeling
 * 
 * Stuff like illegal matrix multiplication and other easily missable developer errors will be throwing exceptions
 */

#ifndef __EXCEPTION_HANDLER__
#define __EXCEPTION_HANDLER__
/**
 * @brief Array for debugger printing.
 */
static const char* exception_name[] = {"Error_none", "Error_not_square_matrix", "Error_illeagal_matrix_multiplication", "Error_illeagal_vector_matrix_multiplication", "Error_vector_and_matrix_size_dont_match", "Error_division_by_zero", "Error_zero_determinant", "Error_determinant_is_nan"};
/**
 * @brief Exception types that will be thrown, enclosed in enum.
 */
enum exception_type {Error_not_square_matrix = 1, Error_illeagal_matrix_multiplication, Error_illeagal_vector_matrix_multiplication, Error_vector_and_matrix_size_dont_match, Error_division_by_zero, Error_zero_determinant, Error_determinant_is_nan};


#endif