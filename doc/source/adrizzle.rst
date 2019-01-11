.. _drizzle:

********************
Image Drizzling Step
********************
The operation of drizzling each input image needs to be performed
twice during processing:

* single drizzle step: this initial step drizzles each image onto
  the final output WCS as separate images

* final drizzle step: this step produces the final combined image
  based on the cosmic-ray masks determined by ``AstroDrizzle``

.. automodule:: drizzlepac.adrizzle
   :members:
   :undoc-members:
