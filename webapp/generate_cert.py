"""
Generate a self-signed TLS certificate/key pair for serving this app over
HTTPS. Pure-Python (uses the `cryptography` package) so it works identically
on Windows/macOS/Linux with no dependency on an external openssl binary.

Usage:
    python generate_cert.py

Writes certs/cert.pem and certs/key.pem (both git-ignored; regenerate per machine).

The certificate is valid for `localhost`/127.0.0.1 PLUS whatever extra
hostnames/IPs are listed in EXTRA_DNS_NAMES / EXTRA_IPS below -- add your
server's LAN hostname/IP there if it needs to be reachable from other
machines on the network (e.g. via an admin-configured direct route, not just
same-machine or an SSH tunnel). Without a matching Subject Alternative Name,
browsers reject the cert outright for that hostname rather than showing the
normal (click-through) self-signed warning.
"""
import datetime
import ipaddress
import pathlib

from cryptography import x509
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.x509.oid import NameOID

CERT_DIR = pathlib.Path(__file__).parent / "certs"

# This deployment is also reachable directly (not just via SSH tunnel) as
# https://bioclaude:8501, per the admin-provided route -- covered here so the
# certificate is valid for that hostname/IP too, not just localhost.
EXTRA_DNS_NAMES = ["bioclaude"]
EXTRA_IPS = ["10.1.100.110"]


def generate_self_signed_cert(cert_path: pathlib.Path, key_path: pathlib.Path, days_valid: int = 825):
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)

    subject = issuer = x509.Name([
        x509.NameAttribute(NameOID.COUNTRY_NAME, "US"),
        x509.NameAttribute(NameOID.ORGANIZATION_NAME, "Local Mass-Transfer Fitting Tool"),
        x509.NameAttribute(NameOID.COMMON_NAME, "localhost"),
    ])

    dns_names = [x509.DNSName("localhost")] + [x509.DNSName(h) for h in EXTRA_DNS_NAMES]
    ip_addrs = [x509.IPAddress(ipaddress.ip_address("127.0.0.1"))] + [
        x509.IPAddress(ipaddress.ip_address(ip)) for ip in EXTRA_IPS
    ]

    now = datetime.datetime.now(datetime.timezone.utc)
    cert = (
        x509.CertificateBuilder()
        .subject_name(subject)
        .issuer_name(issuer)
        .public_key(key.public_key())
        .serial_number(x509.random_serial_number())
        .not_valid_before(now - datetime.timedelta(days=1))
        .not_valid_after(now + datetime.timedelta(days=days_valid))
        .add_extension(
            x509.SubjectAlternativeName(dns_names + ip_addrs),
            critical=False,
        )
        .add_extension(x509.BasicConstraints(ca=True, path_length=None), critical=True)
        .sign(key, hashes.SHA256())
    )

    cert_path.write_bytes(cert.public_bytes(serialization.Encoding.PEM))
    key_path.write_bytes(
        key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.TraditionalOpenSSL,
            encryption_algorithm=serialization.NoEncryption(),
        )
    )


def main():
    CERT_DIR.mkdir(exist_ok=True)
    cert_path = CERT_DIR / "cert.pem"
    key_path = CERT_DIR / "key.pem"

    if cert_path.exists() and key_path.exists():
        print(f"Certificate already exists at {cert_path} -- delete it first to regenerate.")
        return

    generate_self_signed_cert(cert_path, key_path)
    print(f"Self-signed certificate written to: {cert_path}")
    print(f"Private key written to:             {key_path}")
    print(f"Valid for hostnames: localhost, {', '.join(EXTRA_DNS_NAMES)}")
    print(f"Valid for IPs:       127.0.0.1, {', '.join(EXTRA_IPS)}")
    print("\nThese are SELF-SIGNED -- your browser will show a security warning on first")
    print("visit. This is expected for local-only tools; proceed past the warning (e.g.")
    print("'Advanced -> Proceed'). Do not reuse this certificate for anything internet-facing.")


if __name__ == "__main__":
    main()
