# Partially_saturated_groundwater_flow_model_in_julia

![Preview](https://numericalenvironmental.files.wordpress.com/2018/10/theis-curves.jpg)

This julia script is a complete groundwater flow model based on the integral finite difference method. The model is capable of simulating transient groundwater flow in one-, two-, or three dimensions using an arbitrary volume element geometry; flow can be modeled for confined conditions using a fixed saturated hydraulic conductivity and specific storage, per volume element, or unsaturated conditions by solution of the Richards equation (with relative permeability and capillary pressure curves defined as a function of pressure or fluid saturation), or a mixture of both. More background on the model, and some example application problems, are discussed in my blog, https://numericalenvironmental.wordpress.com/2018/10/26/a-variably-saturated-groundwater-flow-model-coded-in-julia-using-the-integral-finite-difference-method/. Note that hydraulic conductivity anisotropy can only be specified with respect to the vertical direction and requires some user input; it is not automatically assigned across connections (an artifact of the integral finite difference method). Also note that unsaturated hydraulic properties are updated only at the end of each time step; there is currently no iterative solution feature for unsaturated flow and large time steps.  

The code is distributed across two modules: jFlow.jl, which contains all of the core numerical model functions and struct definitions, and jFlowIOModule.jl, which handles reading model input files and writing output. Both modules are set up to run under Julia version 1.0. The following text input files are also required; refer to comments in the source code for the definitions/purposes of individual parameters:

* Cells.txt - set up volume element properties and propagate identical cell definitions in one or more dimensions
* CellMod.txt - apply modifications of volume elements after initial definition in Cells.txt (such as when populating properties according to a correlated random field); this step is optional, so comment out the call to the ReadCellMods function within jFlow.jl if it is not going to be used
* gridConnect.txt - set up regular, repeating volume element connections that can be propagated, such as within a grid
* specConnects.txt - set up special volume element connections, such as multiple connections to a single boundary element, wellbore element, etc.
* Sources.txt - specify source terms, organized by volume element index number, volumetric water flux, and stress period
* Params.txt - specify model run parameters, such as simulation end time, starting and maximum time steps, and maximum change in head per time step
* Monitor.txt - provide list of volume element index numbers to be included in a time series output file

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

