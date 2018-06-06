
unsigned int AvailDeviceMemInBytes()
{
   size_t AvailMemoryInBytes;
   size_t TotalMemoryInBytes;

   cudaMemGetInfo(&AvailMemoryInBytes, &TotalMemoryInBytes);

   return AvailMemoryInBytes;
}
