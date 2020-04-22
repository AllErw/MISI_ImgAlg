# MISI_ImgAlg
Imaging and other algorithms for the Multimodal Interventional Sensing and Imaging (MISI) lab.

Version 1.1
Last edited by Erwin Alles on 22 April 2020.



This library contains a collection of image reconstruction algorithms and support functions tailored to all-optical ultrasound imaging, but many can also be applied to other ultrasound or photoacoustic imaging scenarios. The source code is written in C/C++, and is designed to be compiled into a Dynamically Linked Library (dll) to enable its use in MATLAB, LabVIEW and other languages.

Four different image reconstruction algorithms are implemented:

1. Delay & Sum (DAS)
The simplest and most versatile method, which is also highly efficient. This method is widely documented, but one of the more elegant descriptions is found in [https://doi.org/10.1016/j.ultras.2006.07.017]. The particular implementation in this library is commonly referred to as "synthetic aperture focussing" or "dynamic focussing".

2. Delay, Multiply & Sum (DMAS)
This method further exploits the correlation between neighbouring RF time traces to reduce side lobe image artefacts and other image clutter. It does so by not merely summing time-delayed RF signals, but by performing pair-wise multiplications prior to summation, see [https://doi.org/10.1109/TMI.2014.2371235] for details. This method is especially effective in sparse images (point objects, such as point scatterers or needle tips). Note: this method is non-linear in nature, resulting in an image appearance that differs from conventional ultrasound images, and is significantly slower than DAS.

3. Short-Lag Spatial Coherence (SLSC)
This method employs a different method of exploiting coherence between neighbouring RF time traces, where the coherence between the signals transmitted/received by a narrow range ("short lag") of neighbouring elements is explicitly computed and integrated. See [https://doi.org/10.1177/016173461103300203] for details. Like DMAS, SLSC is non-linear, thus resulting in different image appearance, but SLSC tends to result in images that are somewhat closer to conventional ultrasound images in appearance. SLSC is even slower than DMAS, but typically produces more "useful" images.

4. k-space reconstruction
This is one of the most efficient methods implemented. With this method, the RF data is Fourier transformed into a collection of planewaves, which are analytically time-delayed and then inverse transformed to form the image. While initially developed for photoacoustic image reconstruction [https://doi.org/10.1088/0266-5611/23/6/S05], this method also works for pulse-echo ultrasound imaging where a single source and single receiver (that coincide spatially) are translated in periodic increments across a linear (1D) imaging aperture. This scenario is approximately true for 2-fibre needle imaging probes that are mechanically scanned using translation stages.
Notes: 
1) the version implemented is slightly different than the implementation found in k-Wave, in the sense that it (approximately) performs sinc-interpolation rather than any of the methods implemented in MATLAB. 
2) this method requires periodic sampling in both time and space, and is currently only implemented in 2D!

Method 4 is only implemented for the case where a single source and a single receiver spatially coincide, whereas methods 1-3 are implemented for three cases of increasing versatility and complexity:
a) a single source and a single receiver coincide spatially and are traversed across an aperture or arbitrary geometry (1D, 2D or 3D);
b) a single stationary receiver is used to record pulse-echo signals originating from a multitude of sources distributed across an arbitrarily shaped (1D, 2D or 3D) source aperture;
c) an arbitrary number of sources and receivers are used that are positioned fully independently. This is the most versatile method, and can in principle be applied to conventional (i.e., not-all-optical) ultrasound data, but is also somewhat slower than methods a) and b).

All methods assume point-like (both in spatial extent and bandwidth) sources and receivers, and do not take directivities into account.