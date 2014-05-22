.. _drizzle:

**********************************
Step 4 and 8: Drizzling the Images
**********************************
The operation of drizzling each input image needs to be performed
twice during processing:

* single drizzle step: this initial step drizzles each image onto
  the final output WCS as separate images
* final drizzle step: this step produces the final combined image
  based on the cosmic-ray masks determined by `MultiDrizzle`

.. automodule:: drizzlepac.adrizzle
   :members:
   :undoc-members:
