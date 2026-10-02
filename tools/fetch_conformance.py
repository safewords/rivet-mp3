#!/usr/bin/env python3
"""Fetch the MPEG-1/2 audio (Layer I, II, III) conformance streams and their
reference waveforms from ISO's public ISO/IEC 14496-26 (2nd ed.) conformance
package, and unwrap each stream from its MP4 into a plain MPEG audio
elementary stream.

    python tools/fetch_conformance.py DIR

DIR receives `<name>.mpa` (the stream) and `<name>.wav` (ISO's reference
output) for every l1_*, l2_* and l3_* sequence. Layer I and II streams are
taken from the MP4 files (compressedMp4), whose samples are plain MPEG audio
frames; Layer III streams from the native-format files
(compressedMpeg12/add-opt), because the MP4 versions carry Layer III in a
rearranged form rather than as frames. Only the members needed are
read, through HTTP range requests: the reference archive is ~10 GB, the
part used here ~40 MB. Then run

    MP3_CONFORMANCE_DIR=DIR cargo test --release --test conformance -- --nocapture

This script handles data only (ZIP members, MP4 boxes); it contains no audio
decoding.
"""
import io, os, struct, sys, urllib.request, zipfile

BASE = "https://standards.iso.org/iso-iec/14496/-26/ed-2/en/"


class Remote(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos = url, 0
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"))
        self.size = int(head.headers["Content-Length"])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        rng = "bytes=%d-%d" % (self.pos, self.pos + n - 1)
        data = urllib.request.urlopen(urllib.request.Request(self.url, headers={"Range": rng})).read()
        b[: len(data)] = data
        self.pos += len(data)
        return len(data)


def members(archive):
    z = zipfile.ZipFile(io.BufferedReader(Remote(BASE + archive), 1 << 20))
    for info in z.infolist():
        name = os.path.basename(info.filename)
        if info.is_dir() or "_checksum" in info.filename or name[:3] not in ("l1_", "l2_", "l3_"):
            continue
        yield name, lambda info=info: z.read(info)


def boxes(data, start, end):
    while start + 8 <= end:
        size, kind = struct.unpack(">I4s", data[start : start + 8])
        hdr = 8
        if size == 1:
            size = struct.unpack(">Q", data[start + 8 : start + 16])[0]
            hdr = 16
        elif size == 0:
            size = end - start
        yield kind.decode("latin-1"), start + hdr, start + size
        start += size


def find(data, path, start=0, end=None):
    end = len(data) if end is None else end
    for kind, s, e in boxes(data, start, end):
        if kind == path[0]:
            if len(path) == 1:
                return s, e
            skip = {"stsd": 8}.get(kind, 0)
            r = find(data, path[1:], s + skip, e)
            if r:
                return r
    return None


def mp4_samples(data):
    """Concatenate the samples of the first track, in decoding order."""
    stbl = ("moov", "trak", "mdia", "minf", "stbl")
    s, e = find(data, stbl + ("stsz",))
    _, fixed, count = struct.unpack(">III", data[s : s + 12])
    sizes = [fixed] * count if fixed else list(struct.unpack(">%dI" % count, data[s + 12 : s + 12 + 4 * count]))
    r = find(data, stbl + ("stco",))
    if r:
        s, e = r
        n = struct.unpack(">I", data[s + 4 : s + 8])[0]
        offsets = struct.unpack(">%dI" % n, data[s + 8 : s + 8 + 4 * n])
    else:
        s, e = find(data, stbl + ("co64",))
        n = struct.unpack(">I", data[s + 4 : s + 8])[0]
        offsets = struct.unpack(">%dQ" % n, data[s + 8 : s + 8 + 8 * n])
    s, e = find(data, stbl + ("stsc",))
    n = struct.unpack(">I", data[s + 4 : s + 8])[0]
    runs = [struct.unpack(">III", data[s + 8 + 12 * i : s + 20 + 12 * i]) for i in range(n)]
    out, k = bytearray(), 0
    for chunk, off in enumerate(offsets, 1):
        per = [r[1] for r in runs if r[0] <= chunk][-1]
        for _ in range(per):
            if k == len(sizes):
                break
            out += data[off : off + sizes[k]]
            off += sizes[k]
            k += 1
    return bytes(out)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "conformance"
    os.makedirs(out, exist_ok=True)
    for name, read in members("compressedMpeg12.zip"):
        target = os.path.join(out, name[:-4] + ".mpa")
        if not os.path.exists(target):
            open(target, "wb").write(read())
            print("stream", target)
    for name, read in members("compressedMp4.zip"):
        if name.startswith("l3_"):
            continue
        target = os.path.join(out, name[:-4] + ".mpa")
        if not os.path.exists(target):
            open(target, "wb").write(mp4_samples(read()))
            print("stream", target)
    for name, read in members("referencesWav.zip"):
        target = os.path.join(out, name)
        if not os.path.exists(target):
            open(target, "wb").write(read())
            print("reference", target)


if __name__ == "__main__":
    main()
