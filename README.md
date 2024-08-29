Overview
The Gdz library is a C++ template-based library designed for highly efficient manipulation of arbitrarily large integers, with a primary focus on compactness, speed, and flexibility. The library supports various operations, including arithmetic, bitwise manipulation, comparison, construction from various types, and more. It also includes utilities for floating-point conversions, container interfaces, and detailed bit manipulation.  Extensive constexpr implementation to permit compile time Gdz expressions.  Powerful and complete functionality of intialization, conversion, output - all natively supported with zero syntax bloat.  Gdz behaves identically to 1st class built-in types with very few exceptions.  Refactoring to use it will be almost* as simple as doing a global search and replace on 'int'.  Unbelievably most other big integer implementations use strings as the internal numeric representation.  Gdz uses native 32 and 64 bit operations extended to arbitrary precision operating on an array of 32 bit blocks.  Gdz is likely somewhere between 10x to 1000x faster than 'string arithmetic' alternatives - more? 

Features
Arbitrary Precision Integers:

Gdz supports arbitrary precision by templating the number of blocks (words) used to represent the integer. For example, Gdz<4> represents a POD 2's compliment integer using 4 U32 blocks (128 bits) with no additional data members.

Comprehensive Constructor Support:

Supports various construction methods, including default, copy, initializer list, and construction from a single value or multiple homogeneous/heterogeneous values.

Arithmetic and Bitwise Operations:

Full support for arithmetic (addition, subtraction, multiplication, etc.) and bitwise operations (AND, OR, XOR, NOT, shifts, and rotates).

Comparison Operations:

Provides a full suite of comparison operators (==, !=, <, >, <=, >=), supporting comparisons against both Gdz objects and native C++ types.

Floating-Point Conversion:

Gdz can convert large integers to and from floating-point types, with attention to precision and rounding.
Container Interface:

Gdz provides an interface akin to STL containers, with support for iteration, element access, and size queries.
Bit Manipulation Utilities:

Advanced bit manipulation functions, including bit counting, finding least/most significant bits, and reversing bits.
Examples:

//1. Initialization and Construction

#define DO_CREATE
#include "Gdz.hpp"

int main() {

    // Default constructor
    S128 gdz1 = "11111000:01010101:00001111"_g2; //binary, octal, decimal, hex native literals of arbitrary size
    
    // Copy constructor
    S1k gdz2(gdz1);
    
    // Initializer list constructor
    S128 gdz3{0x82345678, 0x23456789, 0x34567890, 0x45678901};
    
    // Single value constructor
    S128 gdz4(1234);
    
    // Multiple homogeneous values constructor
    S128 gdz5(0x12345678, 0x23456789, 0x34567890, 0x45678901);
    
    // Multiple heterogeneous values constructor
    S128 gdz6(0x82345678, (U16)0x9ABC, 0x13579BDF, (U8)0x24);
    

//2. Arithmetic and Bitwise Operations

    Gdz<2> a{1234};
    S256 b{5678};     // Gdz<8> alias
    int c;
    
    // Arithmetic operations
    auto sum = a + b;
    auto diff = a - c;
    auto prod = a * b;
    auto quot = a / b;
    auto rem = c % b;
    auto rem = c + a % b;
    
    // Bitwise operations
    auto and_result = a & b;
    auto or_result = a | c;
    auto xor_result = a ^ b;
    auto not_result = ~a;
    
    // Shifts
    auto left_shift = a << 2;
    auto right_shift = a >> 2;
    
//3. Comparison Operations

    S128 gdz1{1234};
    Gdz<4> gdz2{5678};  //equivalent size
    
    // Comparison operations
    bool eq = gdz1 == gdz2;
    bool neq = gdz1 != gdz2;
    bool lt = gdz1 < gdz2;
    bool gt = gdz1 > 5;
    bool le = gdz1 <= 10.0;
    bool ge = 10 >= gdz2;
    
//4. Floating-Point Conversion

    double f = 1.23456789e10;
    S256 gdz(f);
    
    // Convert back to float
    double f_back = static_cast<double>(gdz);
    
//5. Container Interface

    S128 gdz{0x12345678, 0x23456789, 0x34567890, 0x45678901};
    
    // Iteration using iterators
    for (auto it = gdz.begin(); it != gdz.end(); ++it) {
        std::cout << *it << " ";
    }
    
    // Range-based for loop
    for (auto val : gdz) {
        std::cout << val << " ";
    }
    
    // Element access
    std::cout << gdz[0] << " " << gdz[1] << " ";
    

//6. Flexible and Easy Formatted Output

    S128 binary = "1011011101111011111"_g2;
    std::cout << binary; // output: {5BBDF} (default: hex, 'big endian' ordering)
    
    GdzT.Set( 2 /*radix*/, 0 /*msd,msb ordering*/, 8 /segmenation*/, ":" /*seperator*/, "[\a]" /*wrap*/ );
    std::cout << binary; // output: [101:10111011:11011111] 

    S256 ID = "Hello_World!"_gID; // the '!' is ignored.  any char outside radix range is ignored as inert formatting
    Gdz.Radix( 37 ); // gID token radix
    std::cout << ID; // output: [HEL:LO_WORLD] -  the segmenation is still 8 from prior set
    Gdz.Flags( 3 );  // lsd,lsb ordering
    std::cout << ID; // output: [DLROW_OL:LEH] - any combination of radix, endianess/ordering, segmentation is valid

    //GdzT is the formatting object that simplifies complex configuation yet allows minimal usage syntax
    //multithreading supported wit instances.
    
    return 0;
}

Installation
Simply include "Gdz.hpp" in your project. No additional libraries or dependencies are required. #define DO_CREATE before inclusion into the primary translation unit to instatiate necessary globals.



This README provides a partial guide for understanding and using the Gdz library, including installation, features, and example usage. The library's emphasis is on performance, flexibility, simplicity and elegance making it an ultra-lightweight and completely integrated large integer solution.  In almost all use cases and syntaxes it functions identically to builtin types which is something no other comparable library can claim. Code written with it will be concise, clear and efficient.  The suite of intialization, conversion, output functionality gives it a 'super connected' syntax.  Output of these extremely large numbers is simple yet completely configurable with it being able to create strings in any combination of radix, segmentation and ordering.  
