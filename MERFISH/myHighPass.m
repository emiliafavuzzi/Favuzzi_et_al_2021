function thisimage=myHighPass(image)
thisimage=image-imgaussfilt(image,5);
