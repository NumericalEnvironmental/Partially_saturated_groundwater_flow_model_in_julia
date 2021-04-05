# Partially_saturated_groundwater_flow_model_in_julia

![Preview](https://numericalenvironmental.files.wordpress.com/2018/10/theis-curves.jpg)

April 5, 2021 Update:

jFlow has been modified to (1) work with an externally generated numerical grid/mesh, and (2) accept a linearized relationship between relative permeability and water saturation in the local pressure head on a cell-specific basis, enabling easier simulation of unconfined flow, vertical conduits, and (planned) surface water features. See my blog post [link pending] for more details.

The file set within the repository supports the current example file in my blog. This includes the nodes.csv and connections.csv files, which were created using a cylindrical mesh generator specifically designed to be compatible with jFlow (the mesh generator – a Python script – and its supporting files are included in a subfolder). The “mesh” parameter within the Params.txt input file must be set to “true” to enable reading these files. Otherwise, the default Cells.txt and gridConnections.txt files are used to set up the model grid (older, unrelated setup files for a previous demo problem remain in the folder to illustrate formatting for these).

Separately, the previously unused “b” parameter in the cell struct now serves a purpose: setting the value equal to zero enables the Van Genuchten unsaturated flow properties parameterization for the particular grid cell, whereas setting “b” equal to a value greater than zero corresponds to the cell vertical thickness, which is used in the linearized relative permeability and water saturation pressure head formulation.

Previous documentation:

This julia script is a complete groundwater flow model based on the integral finite difference method. The model is capable of simulating transient groundwater flow in one-, two-, or three dimensions using an arbitrary volume element geometry; flow can be modeled for confined conditions using a fixed saturated hydraulic conductivity and specific storage, per volume element, or unsaturated conditions by solution of the Richards equation (with relative permeability and capillary pressure curves defined as a function of pressure or fluid saturation), or a mixture of both. More background on the model, and some example application problems, are discussed in my blog, https://numericalenvironmental.wordpress.com/2018/10/26/a-variably-saturated-groundwater-flow-model-coded-in-julia-using-the-integral-finite-difference-method/. Note that hydraulic conductivity anisotropy can only be specified with respect to the vertical direction and requires some user input; it is not automatically assigned across connections (an artifact of the integral finite difference method). Also note that unsaturated hydraulic properties are updated only at the end of each time step; there is currently no iterative solution feature for unsaturated flow and large time steps.  

The code is distributed across two modules: jFlow.jl, which contains all of the core numerical model functions and struct definitions, and jFlowIOModule.jl, which handles reading model input files and writing output. Both modules are set up to run under Julia version 1.0. The following text input files are also required; refer to comments in the source code for the definitions/purposes of individual parameters:

* Cells.txt - set up volume element properties and propagate identical cell definitions in one or more dimensions
* CellMod.txt - apply modifications of volume elements after initial definition in Cells.txt (such as when populating properties according to a correlated random field); this step is optional, so comment out the call to the ReadCellMods function within jFlow.jl if it is not going to be used
* gridConnect.txt - set up regular, repeating volume element connections that can be propagated, such as within a grid
* specConnects.txt - set up special volume element connections, such as multiple connections to a single boundary element, wellbore element, etc.
* Sources.txt - specify source terms, organized by volume element index number, volumetric water flux, and stress period
* Drains.txt –indices for those cells with a fluid pressure cap (i.e. drain or seepage face), with corresponding threshold pressures
* Params.txt - specify model run parameters, such as simulation end time, starting and maximum time steps, and maximum change in head per time step
* Monitor.txt - provide list of volume element index numbers to be included in a time series output file

The input files included herein support a vertical two-dimensional unsaturated flow demonstration featuring a seepage face constraint (https://numericalenvironmental.wordpress.com/2019/09/30/an-enhancement-to-a-saturated-unsaturated-groundwater-model-written-in-julia-language/). Note that the CellMod and specGridConnect features are not used for this problem, so the steps that read the corresponding files in the jFlow.jl code are commented out. These input files are included here nonetheless to communicate their required formats, when needed for other applications.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

