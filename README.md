# liton_point_data

Point_data is a simple array class in c++. This class is designed for the storage of point data and implemented with template in c++11.

## Features

### Basic Features
* Support 0, 1, 2 and 3 dimensional point data
* Support components at each point with SOA (structure  of array) layout
* Index with point number which can be negative index
* Support index range and flag check in debug mode
* Provide basic display of array message and data in ascii mode

### Highlight Features
* Allow some virtual point in front of and ahead of the main domain in each dimension
* Provide half point and center point in each dimension
* Easy to get often used range of the domain point
* Support loop and reduce operation with lambda expression

## Usage

### Download and Configuration
You only need to download all folders in `scr`  for using this package.

Before using it in your project, you need to copy these folders into your project directory or somewhere else to make sure that the complier could find it. 

### Compile
This library uses c++11, make sure your compiler support it and turn on the option.

## Contribution
You can just comment in issues or contact me by [email](mailto:luan_ming_yi@126.com) to talk about your idea or demand. Thank you.

## License
[MIT License](https://opensource.org/licenses/MIT)
