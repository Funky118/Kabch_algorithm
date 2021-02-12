\mainpage Introduction to MIPO
## Multipoint Identification of the Position and Orientation of an object

This documentation covers the entire code for the MIPO project.\n\n
When working with manipulators it is often necessary to translate between two coordinate systems when working on an object. The first coordinate system relates to the virtual model saved inside the robot\n
and the second to the real object the robot is somehow manipulating with.\n
This can be solved in a number of ways, but this particular documented code uses the Kabsch algorithm to compute a transformation matrix which when applied to the real objects coordinates will transform them\n
into model coordinates.\n
The word for such a model or a real object matrix is point cloud as it is represented by a matrix of 3d points resembling a cloud.\n
If these two point clouds relate to the same object and if their points match in number and indices, then it is possible to use the Kabsch algorithm to calculate the rotation and translation matrix mapping\n
the model point cloud into the real objects one.\n\n

In the root you will find a py/ folder containing a test script for the algorithm.\n\n

In this documentation you will find all the algorithms used in the computation of the Kabsch algorithm, such as the QR algorithm, Householder reflection or Singular Value Decomposition.\n
I hope you will find this documentation helpful and that the presented code is at least readable enough :)\n
Good luck!

\n\n\n<b>Email:</b> dominik.ricanek@gmail.com
