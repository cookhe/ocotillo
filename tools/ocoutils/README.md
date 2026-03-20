# Install Ocotillo Utilities

Until I set this up as an installable python module, import ocoutils using:

```python
sys.path.append('/path/to/ocotillo/tools/ocoutils')

import ReadOcotillo as oco

# exammple usage to open "wavelengths.bin"
oco.read_waves('ocotillo/output/')
```
