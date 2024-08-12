Overview
The Gdz library is a C++ template-based library designed for highly efficient manipulation of arbitrarily large integers, with a primary focus on compactness, speed, and flexibility. The library supports various operations, including arithmetic, bitwise manipulation, comparison, construction from various types, and more. It also includes utilities for floating-point conversions, container interfaces, and detailed bit manipulation.  Extensive constexpr implementation to permit compile time Gdz expressions.  Powerful and complete functionality of intialization, conversion, output - all natively supported with zero syntax bloat.  Gdz behaves identically to 1st class built-in types with very few exceptions.  Refactoring to use it will be almost* as simple as doing a global search and replace on 'int'.  Try that with any other large integer library ;)

Features
Arbitrary Precision Integers:

Gdz supports arbitrary precision by templating the number of blocks (words) used to represent the integer. For example, Gdz<4> represents a POD 2's compliment integer using 4 U32 blocks with no additional data members.
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
Examples
1. Initialization and Construction
#include "Gdz.hpp"

int main() {
    // Default constructor
    S128 gdz1;
    
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
    
    return 0;
}

2. Arithmetic and Bitwise Operations
#include "Gdz.hpp"

int main() {
    Gdz<2> a{1234};
    Gdz<2> b{5678};
    
    // Arithmetic operations
    auto sum = a + b;
    auto diff = a - b;
    auto prod = a * b;
    auto quot = a / b;
    auto rem = a % b;
    
    // Bitwise operations
    auto and_result = a & b;
    auto or_result = a | b;
    auto xor_result = a ^ b;
    auto not_result = ~a;
    
    // Shifts
    auto left_shift = a << 2;
    auto right_shift = a >> 2;
    
    return 0;
}

3. Comparison Operations
#include "Gdz.hpp"

int main() {
    S128 gdz1{1234};
    Gdz<4> gdz2{5678};  //equivalent 
    
    // Comparison operations
    bool eq = gdz1 == gdz2;
    bool neq = gdz1 != gdz2;
    bool lt = gdz1 < gdz2;
    bool gt = gdz1 > gdz2;
    bool le = gdz1 <= gdz2;
    bool ge = gdz1 >= gdz2;
    
    return 0;
}
4. Floating-Point Conversion
#include "Gdz.hpp"

int main() {
    double f = 1.23456789e10;
    S256 gdz(f);
    
    // Convert back to float
    double f_back = static_cast<double>(gdz);
    
    return 0;
}

5. Container Interface
#include "Gdz.hpp"

int main() {
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
    
    return 0;
}
Installation
Simply include the relevant header files in your project. No additional libraries or dependencies are required. #define DO_CREATE before inclusion into the primary translation unit to instatiate necessary globals.



This README provides a partial guide for understanding and using the Gdz library, including installation, features, and example usage. The library's emphasis is on performance, flexibility, simplicity and elegance making it an ultra-lightweight and completely integrated large integer solution.  In almost all use cases and syntaxes it functions identically to builtin types which is something no other comparable library can claim. Code written with it will be concise, clear and efficient.  The suite of intialization, conversion, output functionality gives it a 'super connected' syntax.  Output of these extremely large numbers is simple yet completely configurable with it being able to create strings in any combination of radix, segmentation and ordering.  
