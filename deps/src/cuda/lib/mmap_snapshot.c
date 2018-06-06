
// This function can create (if nonexisting) and/or
// change mapped file
double *mmap_snapshot(char *snp_path, long long dim)
{
  struct stat snp;
  /* 8 byte for a double precision */
  long snp_size = dim*8;
  /* create an nonexisting output file, if non exist*/
  long fid = open(snp_path,O_RDWR|O_APPEND|O_CREAT,
    S_IRUSR|S_IWUSR|S_IROTH|S_IWOTH);
  /* extend file to reqired file size */
  ftruncate(fid, snp_size);
  /* obtain output file size after extension as an input of memory mapping */
  stat(snp_path, &snp);
  /* create output file memory mapping using mmap */
  double *snp_ptr = (double*) mmap(NULL, snp.st_size,
  PROT_READ|PROT_WRITE, MAP_SHARED, fid, 0);
  close(fid);
  /* check if memory mapping succeeded */
  if(snp_ptr == MAP_FAILED)
  {
    printf("%s\n","snapshot mapping to memory failed");
  }
  return(snp_ptr);
}


void munmap_snapshot(double *snp_ptr, char *snp_path)
{
  struct stat snp;
  stat(snp_path, &snp);
  munmap(snp_ptr, snp.st_size);
}
