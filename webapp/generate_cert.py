"""
Generate a self-signed TLS certificate/key pair for serving this app over
https://localhost. Pure-Python (uses the `cryptography` package) so it works
identically on Windows/macOS/Linux with no dependency on an external openssl binary.

Usage:
    python generate_cert.py

Writes certs/cert.pem and certs/key.pem (both git-ignored; regenerate per machine).
"""
import datetime
import pathlib

from cryptography import x509
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.x509.oid import NameOID

CERT_DIR = pathlib.Path(__file__).parent / "certs"


def generate_self_signed_cert(cert_path: pathlib.Path, key_path: pathlib.Path, days_valid: int = 825):
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)

    subject = issuer = x509.Name([
        x509.NameAttribute(NameOID.COUNTRY_NAME, "US"),
        x509.NameAttribute(NameOID.ORGANIZATION_NAME, "Local Mass-Transfer Fitting Tool"),
        x509.NameAttribute(NameOID.COMMON_NAME, "localhost"),
    ])

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
            x509.SubjectAlternativeName([
                x509.DNSName("localhost"),
                x509.IPAddress(__import__("ipaddress").ip_address("127.0.0.1")),
            ]),
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
    print("\nThese are SELF-SIGNED -- your browser will show a security warning on first")
    print("visit to https://localhost:8501. This is expected for local-only tools; proceed")
    print("past the warning (e.g. 'Advanced -> Proceed to localhost'). Do not reuse this")
    print("certificate for anything internet-facing.")


if __name__ == "__main__":
    main()
