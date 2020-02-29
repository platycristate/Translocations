import numpy as np
import cv2


# taken from AV
def Align(Image1, Image2):
	Image1 = np.array(Image1)
	Image2 = np.array(Image2)
	def Gradient(Image):
		# Calculate the x and y gradients using Sobel operator
		grad_x = cv2.Sobel(Image, cv2.CV_32F, 1, 0, ksize = 3)
		grad_y = cv2.Sobel(Image, cv2.CV_32F, 0, 1, ksize = 3)
		# Combine the two gradients
		grad = cv2.addWeighted(np.absolute(grad_x), 0.5, np.absolute(grad_y), 0.5, 0)
		return grad
		
	Temp1 = cv2.resize(np.float32(Image1), dsize = (0, 0), fx = 1, fy = 1, interpolation = cv2.INTER_LANCZOS4)
	Temp2 = cv2.resize(np.float32(Image2), dsize = (0, 0), fx = 1, fy = 1, interpolation = cv2.INTER_LANCZOS4)
	warp_mode = cv2.MOTION_AFFINE
	WarpMatrix = np.eye(2, 3, dtype = 'float32') 
	criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 5000,  1e-10)
	(cc, WarpMatrix) = cv2.findTransformECC(Gradient(Temp1), Gradient(Temp2), WarpMatrix, warp_mode, criteria, None, 1)
	WarpMatrix[:, 2] = WarpMatrix[:, 2] / 1
	return cv2.warpAffine(Image2, WarpMatrix, Image2.shape[::-1], flags=cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP), WarpMatrix





   	