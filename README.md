This repository *should* allow you to register a model of the spinal cord and vertebrae into OPM sensor space

The way it works is:
1) register a subject mesh (optical scan, MRI, etc) to the OPM space (we have used 3 OPMS as reference points on the left ear, colloarbones and nason)
2) register generic torso model to the subject , which is now in OPM sensor space
3) load in the heart, lungs, spinal cord and vertebra
4) save and use the information!

if you are NOT using one of the provided spinal cord or vetebra models - i recomend having the spinal cord and bone models in the same space as they are for the subject mesh. 
this way you can apply the same transformation applied to the subject, to the cord and bone too. Otherwise you will need to create a custom way of registering in the spinal cord

if you only have a custom spinal cord mesh (already in the same space as the subject mesh) but no bone model, i have provided an example of how to create a bone model too. just remeber to 
apply any transformation you apply to the cord to the bone as well
