# Extended Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

This project utilizes a Kalman filter to estimate the state of a moving object of interest with noisy Lidar and radar measurements. The project obtains RMSE values that are lower than the tolerance outlined in the project rubric (.11 for positional x,y estimations and 0.52 for motion v_x and v_y estimations). 

*This project requires the Udacity Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)*

*This repository includes two files that can be used to set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.*

*Once the install for uWebSocketIO is complete, the main program can be built and run by doing the following from the project top directory.*

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

*Here is the main protocol that main.cpp uses for uWebSocketIO in communicating with the simulator.*

*INPUT: values provided by the simulator to the c++ program
*
*["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)*


*OUTPUT: values provided by the c++ program to the simulator*

*["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]*

---

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

[//]: # (Image References)

[EKF_map]: ./Docs/EKF_map.png "Extended Kalman Filter Algorithm Map"

## Accuracy
The RSME values for px,py,vx,vy vs ground truth settle out to less than 0.11,0.11,0.52 and 0.52 respectively. 

## Program flow
The program follows the flow as outlined in this diagram ![EKF_map][EKF_map]
Where upon the first measurement the position is recorded as reported to initialize the state vector, converting polar coordinates to Cartesian if the first measurement came from Radar, and on subsequent time steps making predictions using the previous values and timestep. After the prediction is made the update step uses different functions depending on the sensor type. Since they see the world differently, the measurements need to be handled differently.

The Lidar measurement using a linear update method whereas the radar uses a non-linear one, necessitating the need for a jacobian matrix which is handled in tools.cpp.

## Code Efficiency
The code strives to be clean and readable, and to avoid loops, and creating unnecessary complexity whether through the use of loops or control flow checks. 

Note: *Italicized* sections and the "Other Important Dependencies" section were written by Udacity. The vast majority of code in main.cpp was written by Udacity.